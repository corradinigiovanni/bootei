// =====================================================================
//  BOOTEI – Bootstrap Ensemble Inference with lexicographic tie-breaking
//
//  Implements χ², Spearman, Mann–Whitney, Kruskal–Wallis using:
//
//    Primary statistic: T†(Z)
//      - one-sided:   T† = T
//      - two-sided:   T† = |T - ref|
//
//    Tie-breaker score: S_B(Z;W) = (1/B) * sum_{b=1}^B T†(Z^{*(b)})
//      (computed with fixed resampling keys W)
//
//    BOOTEI p-value (add-one Monte Carlo):
//      count permutations with
//        T†(πZ) more extreme than T†(Z), OR
//        T†(πZ) tied with T†(Z) AND S_B(πZ; π·W) >= / <= S_B(Z;W)
//      (direction depends on alternative; two-sided uses T† distance)
//
//  Key implementation point (paper-aligned):
//    - Fix resampling keys W once (Sobol/QMC indices or Efron indices).
//    - Reuse the same keys across permutations (coupling).
//    - For sobol_shift, draw one shift once, then keys are fixed.
//
//  Features retained from your engine:
//    - B <= 1: classical permutation p-value (optional midp).
//    - boot_type: "sobol" | "sobol_shift" | "efron"
//    - perm_seed, keep_perm_stats.
//    - For BOOTEI, tie-breaker computed only when needed (ties).
//
//  Compile in R with: Rcpp::sourceCpp("bootei.cpp")
// =====================================================================

#include <Rcpp.h>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <cmath>
#include <functional>
#include <string>
#include <vector>
#include <cstdint>
#include <memory>


using namespace Rcpp;

// -------------------------------------------------------------------
// RNG state guard: restore .Random.seed on exit (no external side effects)
// -------------------------------------------------------------------
class RNGStateGuard {
public:
  explicit RNGStateGuard(bool active_)
    : env(Rcpp::Environment::global_env()),
      active(active_),
      had_seed(false),
      saved_seed(R_NilValue) {

    if (!active) return;

    had_seed = env.exists(".Random.seed");
    if (had_seed) {
      saved_seed = Rf_duplicate(env.get(".Random.seed"));
      R_PreserveObject(saved_seed);
    }
  }

  ~RNGStateGuard() {
    if (!active) return;

    try {
      if (had_seed) {
        env[".Random.seed"] = saved_seed;
      } else {
        if (env.exists(".Random.seed")) env.remove(".Random.seed");
      }
    } catch (...) {
      // never throw from destructor
    }

    if (saved_seed != R_NilValue) R_ReleaseObject(saved_seed);
  }

private:
  Rcpp::Environment env;
  bool active;
  bool had_seed;
  SEXP saved_seed;
};



// -------------------------------------------------------------------
// Bootstrap type
// -------------------------------------------------------------------
enum BootType {
  BOOT_SOBOL = 0,
  BOOT_SOBOL_SHIFT = 1,
  BOOT_EFRON = 2
};

inline BootType parse_boot_type(const std::string &type) {
  if (type == "sobol") return BOOT_SOBOL;
  if (type == "sobol_shift") return BOOT_SOBOL_SHIFT;
  if (type == "efron") return BOOT_EFRON;
  stop("boot_type must be 'sobol', 'sobol_shift', or 'efron'");
  return BOOT_SOBOL;
}

inline std::string describe_boot_type(BootType bt) {
  switch (bt) {
  case BOOT_SOBOL:       return "Sobol QMC (no shift)";
  case BOOT_SOBOL_SHIFT: return "Sobol QMC (random digital shift)";
  case BOOT_EFRON:       return "Efron bootstrap";
  default:               return "Unknown bootstrap";
  }
}

// -------------------------------------------------------------------
// Helpers
// -------------------------------------------------------------------
inline bool is_na_double(double x) { return NumericVector::is_na(x); }

/* average ranks (ties handled) */
NumericVector rank_numeric(const NumericVector &x) {
  int n = x.size();
  std::vector< std::pair<double,int> > ord(n);
  for (int i = 0; i < n; ++i) ord[i] = {x[i], i};
  std::sort(ord.begin(), ord.end(),
            [](const auto& a, const auto& b){ return a.first < b.first; });
  NumericVector rk(n);
  int i = 0;
  while (i < n) {
    int j = i;
    while (j + 1 < n && ord[j+1].first == ord[i].first) ++j;
    double r = (i + j + 2) / 2.0; // 1-based average rank
    for (int k = i; k <= j; ++k) rk[ ord[k].second ] = r;
    i = j + 1;
  }
  return rk;
}

// -------------------------------------------------------------------
// Sobol 1D (base-2), direct i-th point (1-based) + optional shift
// -------------------------------------------------------------------
inline double sobol_1d_raw(uint32_t index) {
  if (index == 0u) return 0.0;
  const int BITS = 32;
  uint32_t X = 0u;
  uint32_t g = index ^ (index >> 1); // Gray code
  for (int k = 0; k < BITS; ++k) {
    if (g & (1u << k)) {
      X ^= (1u << (BITS - 1 - k));
    }
  }
  const double norm = 4294967296.0; // 2^32
  return (double)X / norm;
}

inline double sobol_1d(uint32_t index, double shift) {
  double u = sobol_1d_raw(index);
  u += shift;
  if (u >= 1.0) u -= 1.0;
  return u;
}

// -------------------------------------------------------------------
// Primary transform T† (paper): one-sided => T, two-sided => |T - ref|
// For "less", we keep T† = T and reverse tail direction in comparisons.
// -------------------------------------------------------------------
inline double T_dagger(double T_raw, const std::string &alternative, double ref) {
  if (is_na_double(T_raw)) return NA_REAL;
  if (alternative == "two.sided") return std::fabs(T_raw - ref);
  return T_raw; // greater or less
}

// -------------------------------------------------------------------
// χ² helpers (integer-coded rows/cols)
// -------------------------------------------------------------------
inline double chi2_from_table(const std::vector< std::vector<double> > &tab,
                              const std::vector<double> &rs,
                              const std::vector<double> &cs,
                              double tot) {
  int nr = (int)tab.size();
  int nc = (int)tab[0].size();
  if (tot <= 0.0) return NA_REAL;
  double chi2 = 0.0;
  for (int r = 0; r < nr; ++r) {
    for (int c = 0; c < nc; ++c) {
      double expct = rs[r] * cs[c] / tot;
      if (expct <= 0.0) continue;
      double diff = tab[r][c] - expct;
      chi2 += (diff * diff) / expct;
    }
  }
  return chi2;
}

double chisq_from_pairs_int(const IntegerVector &row_idx,
                            const IntegerVector &col_idx,
                            int nr, int nc) {
  int n = row_idx.size();
  if (n <= 1 || n != col_idx.size()) return NA_REAL;

  std::vector< std::vector<double> > tab(nr, std::vector<double>(nc, 0.0));
  std::vector<double> rs(nr, 0.0), cs(nc, 0.0);
  double tot = 0.0;

  for (int i = 0; i < n; ++i) {
    int r = row_idx[i];
    int c = col_idx[i];
    if (r < 0 || r >= nr || c < 0 || c >= nc) continue;
    tab[r][c] += 1.0;
    rs[r]     += 1.0;
    cs[c]     += 1.0;
    tot       += 1.0;
  }

  return chi2_from_table(tab, rs, cs, tot);
}

// -------------------------------------------------------------------
// Spearman, Mann–Whitney, Kruskal–Wallis
// -------------------------------------------------------------------
double spearman_raw(const NumericVector &x, const NumericVector &y){
  int n = x.size();
  if (n < 3 || n != y.size()) return NA_REAL;

  NumericVector rx = rank_numeric(x);
  NumericVector ry = rank_numeric(y);

  double mx = mean(rx), my = mean(ry);
  double num = 0.0, dx = 0.0, dy = 0.0;
  for (int i = 0; i < n; ++i) {
    double cx = rx[i] - mx;
    double cy = ry[i] - my;
    num += cx * cy;
    dx  += cx * cx;
    dy  += cy * cy;
  }
  if (dx == 0.0 || dy == 0.0) return NA_REAL;
  return num / std::sqrt(dx * dy);
}

double mannwhitney_raw_U(const NumericVector &x, const NumericVector &y){
  int nx = x.size(), ny = y.size();
  if (nx == 0 || ny == 0) return NA_REAL;

  int n = nx + ny;
  NumericVector z(n);
  IntegerVector grp(n);
  for (int i = 0; i < nx; ++i) { z[i] = x[i]; grp[i] = 0; }
  for (int j = 0; j < ny; ++j) { z[nx + j] = y[j]; grp[nx + j] = 1; }

  NumericVector rk = rank_numeric(z);
  double R1 = 0.0;
  for (int i = 0; i < n; ++i) {
    if (grp[i] == 0) R1 += rk[i];
  }
  return R1 - (double)nx * (nx + 1) / 2.0;
}

std::function<double(const NumericVector&, const IntegerVector&)>
  kruskal_factory(const IntegerVector &g_ref) {
    std::unordered_set<int> seen;
    std::vector<int> levels;
    levels.reserve(g_ref.size());

    for (int i = 0; i < g_ref.size(); ++i) {
      if (IntegerVector::is_na(g_ref[i])) continue;
      int v = g_ref[i];
      if (!seen.count(v)) { seen.insert(v); levels.push_back(v); }
    }
    std::sort(levels.begin(), levels.end());

    int K = (int)levels.size();
    if (K < 2) {
      return [](const NumericVector&, const IntegerVector&){ return NA_REAL; };
    }

    std::unordered_map<int,int> gmap;
    for (int j = 0; j < K; ++j) gmap[ levels[j] ] = j;

    return [=](const NumericVector &x, const IntegerVector &g){
      int nloc = x.size();
      if (nloc != g.size()) return NA_REAL;

      std::vector<double> x2;
      std::vector<int> grp2;
      x2.reserve(nloc);
      grp2.reserve(nloc);

      for (int i = 0; i < nloc; ++i) {
        if (NumericVector::is_na(x[i]) || IntegerVector::is_na(g[i])) continue;
        auto it = gmap.find(g[i]);
        if (it == gmap.end()) continue;
        x2.push_back(x[i]);
        grp2.push_back(it->second);
      }

      int N = (int)x2.size();
      if (N < 2) return NA_REAL;

      NumericVector xvec(N);
      for (int i = 0; i < N; ++i) xvec[i] = x2[i];
      NumericVector rk = rank_numeric(xvec);

      std::vector<double> Rsum(K, 0.0);
      std::vector<int> nj(K, 0);
      for (int i = 0; i < N; ++i) {
        int gidx = grp2[i];
        Rsum[gidx] += rk[i];
        nj[gidx]   += 1;
      }

      int k_eff = 0;
      for (int j = 0; j < K; ++j) if (nj[j] > 0) k_eff++;
      if (k_eff < 2) return NA_REAL;

      double Hsum = 0.0;
      for (int j = 0; j < K; ++j) {
        if (nj[j] == 0) continue;
        Hsum += (Rsum[j] * Rsum[j]) / (double)nj[j];
      }

      double Nd = (double)N;
      double H = (12.0 / (Nd * (Nd + 1.0))) * Hsum - 3.0 * (Nd + 1.0);
      return H;
    };
  }

// -------------------------------------------------------------------
// Resampling keys: generate once and reuse across permutations (coupling)
// -------------------------------------------------------------------
std::vector<int> make_keys_paired(int n, int B, BootType btype, double shift) {
  std::vector<int> idx;
  idx.resize((size_t)B * (size_t)n);

  if (B <= 1) return idx;

  if (btype == BOOT_EFRON) {
    RNGScope scope;
    for (int t = 0; t < B*n; ++t) {
      double u = R::runif(0.0, 1.0);
      int j = (int)std::floor(n * u);
      if (j < 0) j = 0;
      if (j >= n) j = n - 1;
      idx[(size_t)t] = j;
    }
  } else {
    // Sobol (shift may be 0 or random but fixed)
    for (int t = 0; t < B*n; ++t) {
      uint32_t ind = (uint32_t)((t + 1) & 0xFFFFFFFFu); // 1-based
      double u = sobol_1d(ind, shift);
      int j = (int)std::floor(n * u);
      if (j < 0) j = 0;
      if (j >= n) j = n - 1;
      idx[(size_t)t] = j;
    }
  }
  return idx;
}

struct UnpairedKeys {
  std::vector<int> idx_x;
  std::vector<int> idx_y;
};

UnpairedKeys make_keys_unpaired(int nx, int ny, int B, BootType btype,
                                double shift_x, double shift_y) {
  UnpairedKeys K;
  K.idx_x.resize((size_t)B * (size_t)nx);
  K.idx_y.resize((size_t)B * (size_t)ny);

  if (B <= 1) return K;

  if (btype == BOOT_EFRON) {
    RNGScope scope;
    for (int t = 0; t < B*nx; ++t) {
      double u = R::runif(0.0, 1.0);
      int j = (int)std::floor(nx * u);
      if (j < 0) j = 0;
      if (j >= nx) j = nx - 1;
      K.idx_x[(size_t)t] = j;
    }
    for (int t = 0; t < B*ny; ++t) {
      double u = R::runif(0.0, 1.0);
      int j = (int)std::floor(ny * u);
      if (j < 0) j = 0;
      if (j >= ny) j = ny - 1;
      K.idx_y[(size_t)t] = j;
    }
  } else {
    const uint32_t OFFSET = 123456789u;
    for (int t = 0; t < B*nx; ++t) {
      uint32_t ind = (uint32_t)((t + 1) & 0xFFFFFFFFu);
      double u = sobol_1d(ind, shift_x);
      int j = (int)std::floor(nx * u);
      if (j < 0) j = 0;
      if (j >= nx) j = nx - 1;
      K.idx_x[(size_t)t] = j;
    }
    for (int t = 0; t < B*ny; ++t) {
      uint32_t ind = (uint32_t)((t + 1) & 0xFFFFFFFFu);
      double u = sobol_1d(ind + OFFSET, shift_y);
      int j = (int)std::floor(ny * u);
      if (j < 0) j = 0;
      if (j >= ny) j = ny - 1;
      K.idx_y[(size_t)t] = j;
    }
  }
  return K;
}

// -------------------------------------------------------------------
// Bagged tie-breaker S_B: average of T† over bootstrap resamples
// keys are fixed (coupled) and reused across permutations.
// -------------------------------------------------------------------
template<typename VECX, typename VECY, typename STATFN>
double bagged_paired_Tdag(const VECX &x, const VECY &y,
                          STATFN raw_statfn,
                          const std::vector<int> &keys, int B,
                          const std::string &alternative, double ref) {
  int n = x.size();
  if (n != (int)y.size() || n <= 1) return NA_REAL;
  if (B <= 1) {
    double Traw = raw_statfn(x, y);
    return T_dagger(Traw, alternative, ref);
  }

  VECX xb(n);
  VECY yb(n);

  double acc = 0.0;
  int ok = 0;

  for (int b = 0; b < B; ++b) {
    const int base = b * n;
    for (int s = 0; s < n; ++s) {
      int idx = keys[(size_t)base + (size_t)s];
      xb[s] = x[idx];
      yb[s] = y[idx];
    }
    double Traw = raw_statfn(xb, yb);
    double Td   = T_dagger(Traw, alternative, ref);
    if (!is_na_double(Td)) { acc += Td; ++ok; }
  }

  if (ok == 0) {
    double Traw = raw_statfn(x, y);
    return T_dagger(Traw, alternative, ref);
  }
  return acc / (double)ok;
}

// Specialised bagging for χ² on integer pairs (avoid building resample vectors)
double bagged_chisq_Tdag_int(const IntegerVector &row_idx,
                             const IntegerVector &col_idx,
                             int nr, int nc,
                             const std::vector<int> &keys, int B) {
  int n = row_idx.size();
  if (n != col_idx.size() || n <= 1) return NA_REAL;
  if (B <= 1) return chisq_from_pairs_int(row_idx, col_idx, nr, nc); // already one-sided

  std::vector< std::vector<double> > tab(nr, std::vector<double>(nc, 0.0));
  std::vector<double> rs(nr, 0.0), cs(nc, 0.0);

  double acc = 0.0;
  int ok = 0;

  for (int b = 0; b < B; ++b) {
    for (int r = 0; r < nr; ++r) {
      std::fill(tab[r].begin(), tab[r].end(), 0.0);
      rs[r] = 0.0;
    }
    std::fill(cs.begin(), cs.end(), 0.0);
    double tot = 0.0;

    const int base = b * n;
    for (int s = 0; s < n; ++s) {
      int idx = keys[(size_t)base + (size_t)s];
      int r = row_idx[idx];
      int c = col_idx[idx];
      if (r < 0 || r >= nr || c < 0 || c >= nc) continue;
      tab[r][c] += 1.0;
      rs[r]     += 1.0;
      cs[c]     += 1.0;
      tot       += 1.0;
    }

    double Td = chi2_from_table(tab, rs, cs, tot);
    if (!is_na_double(Td)) { acc += Td; ++ok; }
  }

  if (ok == 0) return chisq_from_pairs_int(row_idx, col_idx, nr, nc);
  return acc / (double)ok;
}

// Unpaired bagging (Mann–Whitney)
double bagged_unpaired_Tdag(const NumericVector &x,
                            const NumericVector &y,
                            const std::function<double(const NumericVector&, const NumericVector&)> &raw_statfn,
                            const UnpairedKeys &K, int B,
                            const std::string &alternative, double ref) {
  int nx = x.size(), ny = y.size();
  if (nx <= 0 || ny <= 0) return NA_REAL;
  if (B <= 1) {
    double Traw = raw_statfn(x, y);
    return T_dagger(Traw, alternative, ref);
  }

  NumericVector xb(nx), yb(ny);
  double acc = 0.0;
  int ok = 0;

  for (int b = 0; b < B; ++b) {
    const int base_x = b * nx;
    const int base_y = b * ny;

    for (int s = 0; s < nx; ++s) xb[s] = x[ K.idx_x[(size_t)base_x + (size_t)s] ];
    for (int s = 0; s < ny; ++s) yb[s] = y[ K.idx_y[(size_t)base_y + (size_t)s] ];

    double Traw = raw_statfn(xb, yb);
    double Td   = T_dagger(Traw, alternative, ref);
    if (!is_na_double(Td)) { acc += Td; ++ok; }
  }

  if (ok == 0) {
    double Traw = raw_statfn(x, y);
    return T_dagger(Traw, alternative, ref);
  }
  return acc / (double)ok;
}

// -------------------------------------------------------------------
// Classical permutation p-value (B <= 1): optional midp
// -------------------------------------------------------------------
template<typename VECX, typename VECY, typename STATFN>
double perm_classical_paired(const VECX &x, const VECY &y, STATFN raw_statfn,
                             double T0_dag,
                             int R,
                             const std::string &alternative,
                             double ref,
                             bool midp,
                             bool keep_perm_stats,
                             NumericVector &perm_primary) {
  if (R <= 0 || is_na_double(T0_dag)) {
    perm_primary = NumericVector(0);
    return NA_REAL;
  }

  RNGScope scope;
  int n = x.size();

  int ge = 0, gt = 0, eq = 0, valid = 0;
  const double eps = 1e-12;

  if (keep_perm_stats) perm_primary = NumericVector(R, NumericVector::get_na());
  else perm_primary = NumericVector(0);

  for (int r = 0; r < R; ++r) {
    IntegerVector perm = Rcpp::sample(n, n, false);
    VECY yp(n);
    for (int i = 0; i < n; ++i) yp[i] = y[perm[i] - 1];

    double Traw = raw_statfn(x, yp);
    double Td   = T_dagger(Traw, alternative, ref);
    if (keep_perm_stats) perm_primary[r] = Td;

    if (is_na_double(Td)) continue;
    ++valid;

    if (alternative == "greater") {
      if (Td >= T0_dag - eps) ge++;
      if (Td >  T0_dag + eps) gt++;
      else if (std::fabs(Td - T0_dag) <= eps) eq++;
    } else if (alternative == "less") {
      if (Td <= T0_dag + eps) ge++;
      if (Td <  T0_dag - eps) gt++;
      else if (std::fabs(Td - T0_dag) <= eps) eq++;
    } else { // two.sided: Td already |T- ref|, use greater tail
      if (Td >= T0_dag - eps) ge++;
      if (Td >  T0_dag + eps) gt++;
      else if (std::fabs(Td - T0_dag) <= eps) eq++;
    }
  }

  if (valid == 0) return NA_REAL;
  if (midp) return ((double)gt + 0.5 * (double)eq) / (double)valid;
  return (double)(ge + 1) / (double)(valid + 1);
}

double perm_classical_unpaired(const NumericVector &x, const NumericVector &y,
                               const std::function<double(const NumericVector&, const NumericVector&)> &raw_statfn,
                               double T0_dag,
                               int R,
                               const std::string &alternative,
                               double ref,
                               bool midp,
                               bool keep_perm_stats,
                               NumericVector &perm_primary) {
  if (R <= 0 || is_na_double(T0_dag)) {
    perm_primary = NumericVector(0);
    return NA_REAL;
  }

  RNGScope scope;
  int nx = x.size(), ny = y.size();
  int N  = nx + ny;
  if (N == 0) {
    perm_primary = NumericVector(0);
    return NA_REAL;
  }

  NumericVector pool(N);
  for (int i = 0; i < nx; ++i) pool[i]    = x[i];
  for (int j = 0; j < ny; ++j) pool[nx+j] = y[j];

  int ge = 0, gt = 0, eq = 0, valid = 0;
  const double eps = 1e-12;

  if (keep_perm_stats) perm_primary = NumericVector(R, NumericVector::get_na());
  else perm_primary = NumericVector(0);

  for (int r = 0; r < R; ++r) {
    IntegerVector perm = Rcpp::sample(N, N, false);
    NumericVector xg(nx), yg(ny);
    for (int i = 0; i < nx; ++i) xg[i] = pool[perm[i] - 1];
    for (int j = 0; j < ny; ++j) yg[j] = pool[perm[nx + j] - 1];

    double Traw = raw_statfn(xg, yg);
    double Td   = T_dagger(Traw, alternative, ref);
    if (keep_perm_stats) perm_primary[r] = Td;

    if (is_na_double(Td)) continue;
    ++valid;

    if (alternative == "greater") {
      if (Td >= T0_dag - eps) ge++;
      if (Td >  T0_dag + eps) gt++;
      else if (std::fabs(Td - T0_dag) <= eps) eq++;
    } else if (alternative == "less") {
      if (Td <= T0_dag + eps) ge++;
      if (Td <  T0_dag - eps) gt++;
      else if (std::fabs(Td - T0_dag) <= eps) eq++;
    } else {
      if (Td >= T0_dag - eps) ge++;
      if (Td >  T0_dag + eps) gt++;
      else if (std::fabs(Td - T0_dag) <= eps) eq++;
    }
  }

  if (valid == 0) return NA_REAL;
  if (midp) return ((double)gt + 0.5 * (double)eq) / (double)valid;
  return (double)(ge + 1) / (double)(valid + 1);
}

// χ² classical permutation (greater-tailed only)
double perm_classical_chisq_int(const IntegerVector &row_idx,
                                const IntegerVector &col_idx,
                                int nr, int nc,
                                double T0,
                                int R,
                                bool midp,
                                bool keep_perm_stats,
                                NumericVector &perm_primary) {
  if (R <= 0 || is_na_double(T0)) {
    perm_primary = NumericVector(0);
    return NA_REAL;
  }

  RNGScope scope;
  int n = row_idx.size();
  if (n != col_idx.size()) {
    perm_primary = NumericVector(0);
    return NA_REAL;
  }

  int ge = 0, gt = 0, eq = 0, valid = 0;
  const double eps = 1e-12;

  IntegerVector col_perm(n);

  if (keep_perm_stats) perm_primary = NumericVector(R, NumericVector::get_na());
  else perm_primary = NumericVector(0);

  for (int r = 0; r < R; ++r) {
    IntegerVector perm = Rcpp::sample(n, n, false);
    for (int i = 0; i < n; ++i) col_perm[i] = col_idx[perm[i] - 1];

    double Td = chisq_from_pairs_int(row_idx, col_perm, nr, nc);
    if (keep_perm_stats) perm_primary[r] = Td;

    if (is_na_double(Td)) continue;
    ++valid;

    if (Td >= T0 - eps) ge++;
    if (Td >  T0 + eps) gt++;
    else if (std::fabs(Td - T0) <= eps) eq++;
  }

  if (valid == 0) return NA_REAL;
  if (midp) return ((double)gt + 0.5 * (double)eq) / (double)valid;
  return (double)(ge + 1) / (double)(valid + 1);
}

// -------------------------------------------------------------------
// BOOTEI permutation p-value (B > 1): lexicographic tie-breaking
// -------------------------------------------------------------------
template<typename VECX, typename VECY, typename STATFN>
double perm_bootei_paired(const VECX &x, const VECY &y,
                          STATFN raw_statfn,
                          double T0_dag, double S0,
                          const std::vector<int> &keys, int B,
                          int R,
                          const std::string &alternative,
                          double ref,
                          bool keep_perm_stats,
                          NumericVector &perm_primary,
                          NumericVector &perm_tie) {
  if (R <= 0 || is_na_double(T0_dag) || is_na_double(S0)) {
    perm_primary = NumericVector(0);
    perm_tie     = NumericVector(0);
    return NA_REAL;
  }

  RNGScope scope;
  int n = x.size();

  int tail = 0, valid = 0;
  const double eps = 1e-12;

  if (keep_perm_stats) {
    perm_primary = NumericVector(R, NumericVector::get_na());
    perm_tie     = NumericVector(R, NumericVector::get_na());
  } else {
    perm_primary = NumericVector(0);
    perm_tie     = NumericVector(0);
  }

  for (int r = 0; r < R; ++r) {
    IntegerVector perm = Rcpp::sample(n, n, false);
    VECY yp(n);
    for (int i = 0; i < n; ++i) yp[i] = y[perm[i] - 1];

    double Traw = raw_statfn(x, yp);
    double Td   = T_dagger(Traw, alternative, ref);
    if (keep_perm_stats) perm_primary[r] = Td;

    if (is_na_double(Td)) continue;
    ++valid;

    bool is_more_extreme = false;
    bool is_tied = (std::fabs(Td - T0_dag) <= eps);

    if (alternative == "greater") {
      if (Td > T0_dag + eps) is_more_extreme = true;
    } else if (alternative == "less") {
      if (Td < T0_dag - eps) is_more_extreme = true;
    } else { // two.sided: Td already distance, greater tail
      if (Td > T0_dag + eps) is_more_extreme = true;
    }

    if (is_more_extreme) {
      tail++;
      continue;
    }

    if (!is_tied) continue;

    // tie: compute S only now
    double Sr = bagged_paired_Tdag(x, yp, raw_statfn, keys, B, alternative, ref);
    if (keep_perm_stats) perm_tie[r] = Sr;
    if (is_na_double(Sr)) continue;

    if (alternative == "greater" || alternative == "two.sided") {
      if (Sr >= S0 - eps) tail++;
    } else { // less
      if (Sr <= S0 + eps) tail++;
    }
  }

  if (valid == 0) return NA_REAL;
  return (double)(tail + 1) / (double)(valid + 1);
}

double perm_bootei_unpaired(const NumericVector &x, const NumericVector &y,
                            const std::function<double(const NumericVector&, const NumericVector&)> &raw_statfn,
                            double T0_dag, double S0,
                            const UnpairedKeys &K, int B,
                            int R,
                            const std::string &alternative,
                            double ref,
                            bool keep_perm_stats,
                            NumericVector &perm_primary,
                            NumericVector &perm_tie) {
  if (R <= 0 || is_na_double(T0_dag) || is_na_double(S0)) {
    perm_primary = NumericVector(0);
    perm_tie     = NumericVector(0);
    return NA_REAL;
  }

  RNGScope scope;
  int nx = x.size(), ny = y.size();
  int N  = nx + ny;
  if (N == 0) {
    perm_primary = NumericVector(0);
    perm_tie     = NumericVector(0);
    return NA_REAL;
  }

  NumericVector pool(N);
  for (int i = 0; i < nx; ++i) pool[i] = x[i];
  for (int j = 0; j < ny; ++j) pool[nx + j] = y[j];

  int tail = 0, valid = 0;
  const double eps = 1e-12;

  if (keep_perm_stats) {
    perm_primary = NumericVector(R, NumericVector::get_na());
    perm_tie     = NumericVector(R, NumericVector::get_na());
  } else {
    perm_primary = NumericVector(0);
    perm_tie     = NumericVector(0);
  }

  for (int r = 0; r < R; ++r) {
    IntegerVector perm = Rcpp::sample(N, N, false);
    NumericVector xg(nx), yg(ny);
    for (int i = 0; i < nx; ++i) xg[i] = pool[perm[i] - 1];
    for (int j = 0; j < ny; ++j) yg[j] = pool[perm[nx + j] - 1];

    double Traw = raw_statfn(xg, yg);
    double Td   = T_dagger(Traw, alternative, ref);
    if (keep_perm_stats) perm_primary[r] = Td;

    if (is_na_double(Td)) continue;
    ++valid;

    bool is_more_extreme = false;
    bool is_tied = (std::fabs(Td - T0_dag) <= eps);

    if (alternative == "greater") {
      if (Td > T0_dag + eps) is_more_extreme = true;
    } else if (alternative == "less") {
      if (Td < T0_dag - eps) is_more_extreme = true;
    } else {
      if (Td > T0_dag + eps) is_more_extreme = true;
    }

    if (is_more_extreme) {
      tail++;
      continue;
    }
    if (!is_tied) continue;

    double Sr = bagged_unpaired_Tdag(xg, yg, raw_statfn, K, B, alternative, ref);
    if (keep_perm_stats) perm_tie[r] = Sr;
    if (is_na_double(Sr)) continue;

    if (alternative == "greater" || alternative == "two.sided") {
      if (Sr >= S0 - eps) tail++;
    } else {
      if (Sr <= S0 + eps) tail++;
    }
  }

  if (valid == 0) return NA_REAL;
  return (double)(tail + 1) / (double)(valid + 1);
}

double perm_bootei_chisq_int(const IntegerVector &row_idx,
                             const IntegerVector &col_idx,
                             int nr, int nc,
                             double T0, double S0,
                             const std::vector<int> &keys, int B,
                             int R,
                             bool keep_perm_stats,
                             NumericVector &perm_primary,
                             NumericVector &perm_tie) {
  if (R <= 0 || is_na_double(T0) || is_na_double(S0)) {
    perm_primary = NumericVector(0);
    perm_tie     = NumericVector(0);
    return NA_REAL;
  }

  RNGScope scope;
  int n = row_idx.size();
  if (n != col_idx.size()) {
    perm_primary = NumericVector(0);
    perm_tie     = NumericVector(0);
    return NA_REAL;
  }

  int tail = 0, valid = 0;
  const double eps = 1e-12;
  IntegerVector col_perm(n);

  if (keep_perm_stats) {
    perm_primary = NumericVector(R, NumericVector::get_na());
    perm_tie     = NumericVector(R, NumericVector::get_na());
  } else {
    perm_primary = NumericVector(0);
    perm_tie     = NumericVector(0);
  }

  for (int r = 0; r < R; ++r) {
    IntegerVector perm = Rcpp::sample(n, n, false);
    for (int i = 0; i < n; ++i) col_perm[i] = col_idx[perm[i] - 1];

    double Td = chisq_from_pairs_int(row_idx, col_perm, nr, nc);
    if (keep_perm_stats) perm_primary[r] = Td;

    if (is_na_double(Td)) continue;
    ++valid;

    if (Td > T0 + eps) {
      tail++;
      continue;
    }

    if (std::fabs(Td - T0) > eps) continue;

    double Sr = bagged_chisq_Tdag_int(row_idx, col_perm, nr, nc, keys, B);
    if (keep_perm_stats) perm_tie[r] = Sr;
    if (is_na_double(Sr)) continue;

    if (Sr >= S0 - eps) tail++;
  }

  if (valid == 0) return NA_REAL;
  return (double)(tail + 1) / (double)(valid + 1);
}

// -------------------------------------------------------------------
// Exported wrapper
// -------------------------------------------------------------------

// [[Rcpp::export]]
List bootei_cpp(SEXP x, SEXP y,
            std::string test = "chisq",
            int B = 100,
            int R = 1000,
            std::string alternative = "two.sided",
            double perm_seed = NA_REAL,
            bool midp = false,                 // only used when B <= 1
            std::string boot_type = "sobol",
            bool keep_perm_stats = false) {

  if (alternative != "two.sided" &&
      alternative != "greater" &&
      alternative != "less") {
    stop("alternative must be 'two.sided', 'greater', or 'less'");
  }

  BootType btype = parse_boot_type(boot_type);

  // Set seed locally (no effect on caller RNG state)
  // Set seed locally (no effect on caller RNG state)
  const bool do_seed = !Rcpp::NumericVector::is_na(perm_seed);
  RNGStateGuard rng_guard(do_seed);

  if (do_seed) {
    Environment base = Environment::namespace_env("base");
    Function set_seed = base["set.seed"];
    set_seed(perm_seed);
  }



  // ---------------------------------------------------------------
  // χ² test (independence)
  // ---------------------------------------------------------------
  if (test == "chisq") {
    CharacterVector xc(x), yc(y);
    if (xc.size() != yc.size()) stop("For chisq test x and y must have equal length.");
    int n = xc.size();
    if (n <= 1) {
      return List::create(_["statistic_raw"]=NA_REAL, _["statistic"]=NA_REAL,
                          _["tie_breaker"]=NA_REAL, _["p.value"]=NA_REAL,
                          _["method"]="BOOTEI χ² test", _["alternative"]="greater");
    }

    // map levels -> 0..nr-1 / 0..nc-1
    std::unordered_map<std::string,int> rmap, cmap;
    rmap.reserve(n); cmap.reserve(n);
    for (int i = 0; i < n; ++i) {
      if (CharacterVector::is_na(xc[i]) || CharacterVector::is_na(yc[i])) continue;
      rmap[ as<std::string>(xc[i]) ];
      cmap[ as<std::string>(yc[i]) ];
    }
    int id = 0; for (auto &kv : rmap) kv.second = id++; int nr = id;
    id = 0; for (auto &kv : cmap) kv.second = id++; int nc = id;

    if (nr < 2 || nc < 2) {
      return List::create(_["statistic_raw"]=NA_REAL, _["statistic"]=NA_REAL,
                          _["tie_breaker"]=NA_REAL, _["p.value"]=NA_REAL,
                          _["method"]="BOOTEI χ² test", _["alternative"]="greater");
    }

    IntegerVector row_idx(n), col_idx(n);
    for (int i = 0; i < n; ++i) {
      if (CharacterVector::is_na(xc[i]) || CharacterVector::is_na(yc[i])) {
        row_idx[i] = -1; col_idx[i] = -1;
      } else {
        row_idx[i] = rmap[ as<std::string>(xc[i]) ];
        col_idx[i] = cmap[ as<std::string>(yc[i]) ];
      }
    }

    // shifts for sobol_shift (fixed once)
    double shift = 0.0;
    if (btype == BOOT_SOBOL_SHIFT) {
      RNGScope scope;
      shift = R::runif(0.0, 1.0);
    }

    // observed raw + primary
    double Traw0 = chisq_from_pairs_int(row_idx, col_idx, nr, nc);
    double T0    = Traw0; // greater-tailed, so T†=T

    // classical vs BOOTEI
    NumericVector perm_primary, perm_tie;

    if (B <= 1) {
      double pval = perm_classical_chisq_int(row_idx, col_idx, nr, nc, T0, R, midp,
                                             keep_perm_stats, perm_primary);

      List out = List::create(
        _["statistic_raw"] = Traw0,
        _["statistic"]     = T0,
        _["tie_breaker"]   = NA_REAL,
        _["p.value"]       = pval,
        _["method"]        = "Permutation χ² test (B=1)",
        _["alternative"]   = "greater"
      );
      if (keep_perm_stats && perm_primary.size() > 0) out["perm_primary"] = perm_primary;
      return out;
    }

    // BOOTEI: fixed keys and tie-breaker
    std::vector<int> keys = make_keys_paired(n, B, btype, shift);
    double S0 = bagged_chisq_Tdag_int(row_idx, col_idx, nr, nc, keys, B);

    double pval = perm_bootei_chisq_int(row_idx, col_idx, nr, nc,
                                        T0, S0, keys, B, R,
                                        keep_perm_stats, perm_primary, perm_tie);

    std::string method_str =
      "BOOTEI χ² test (tie-breaker via " + describe_boot_type(btype) + ")";

    List out = List::create(
      _["statistic_raw"] = Traw0,
      _["statistic"]     = T0,
      _["tie_breaker"]   = S0,
      _["p.value"]       = pval,
      _["method"]        = method_str,
      _["alternative"]   = "greater"
    );
    if (keep_perm_stats) {
      if (perm_primary.size() > 0) out["perm_primary"] = perm_primary;
      if (perm_tie.size() > 0)     out["perm_tie"]     = perm_tie;
    }
    return out;
  }

  // ---------------------------------------------------------------
  // Spearman
  // ---------------------------------------------------------------
  if (test == "spearman") {
    NumericVector xn(x), yn(y);
    if (xn.size() != yn.size()) stop("For spearman test x and y must have equal length.");

    auto raw_statfn = [](const NumericVector &a, const NumericVector &b){
      return spearman_raw(a, b);
    };

    // ref for two-sided: 0
    const double ref = 0.0;

    // shifts
    double shift = 0.0;
    if (btype == BOOT_SOBOL_SHIFT) {
      RNGScope scope;
      shift = R::runif(0.0, 1.0);
    }

    double Traw0 = raw_statfn(xn, yn);
    double T0    = T_dagger(Traw0, alternative, ref);

    NumericVector perm_primary, perm_tie;

    if (B <= 1) {
      double pval = perm_classical_paired(xn, yn, raw_statfn, T0, R, alternative, ref,
                                          midp, keep_perm_stats, perm_primary);

      List out = List::create(
        _["statistic_raw"] = Traw0,
        _["statistic"]     = T0,
        _["tie_breaker"]   = NA_REAL,
        _["p.value"]       = pval,
        _["method"]        = "Permutation Spearman test (B=1)",
        _["alternative"]   = alternative
      );
      if (keep_perm_stats && perm_primary.size() > 0) out["perm_primary"] = perm_primary;
      return out;
    }

    std::vector<int> keys = make_keys_paired(xn.size(), B, btype, shift);
    double S0 = bagged_paired_Tdag(xn, yn, raw_statfn, keys, B, alternative, ref);

    double pval = perm_bootei_paired(xn, yn, raw_statfn,
                                     T0, S0, keys, B, R,
                                     alternative, ref,
                                     keep_perm_stats, perm_primary, perm_tie);

    std::string method_str =
      "BOOTEI Spearman test (tie-breaker via " + describe_boot_type(btype) + ")";

    List out = List::create(
      _["statistic_raw"] = Traw0,
      _["statistic"]     = T0,
      _["tie_breaker"]   = S0,
      _["p.value"]       = pval,
      _["method"]        = method_str,
      _["alternative"]   = alternative
    );
    if (keep_perm_stats) {
      if (perm_primary.size() > 0) out["perm_primary"] = perm_primary;
      if (perm_tie.size() > 0)     out["perm_tie"]     = perm_tie;
    }
    return out;
  }

  // ---------------------------------------------------------------
  // Mann–Whitney (unpaired)
  // ---------------------------------------------------------------
  if (test == "mannwhitney") {
    NumericVector x1(x), y1(y);

    auto raw_statfn = [](const NumericVector &a, const NumericVector &b){
      return mannwhitney_raw_U(a, b);
    };

    const double ref = (double)x1.size() * (double)y1.size() / 2.0;

    double shift_x = 0.0, shift_y = 0.0;
    if (btype == BOOT_SOBOL_SHIFT) {
      RNGScope scope;
      shift_x = R::runif(0.0, 1.0);
      shift_y = R::runif(0.0, 1.0);
    }

    double Traw0 = raw_statfn(x1, y1);
    double T0    = T_dagger(Traw0, alternative, ref);

    NumericVector perm_primary, perm_tie;

    if (B <= 1) {
      double pval = perm_classical_unpaired(x1, y1, raw_statfn, T0, R, alternative, ref,
                                            midp, keep_perm_stats, perm_primary);

      List out = List::create(
        _["statistic_raw"] = Traw0,
        _["statistic"]     = T0,
        _["tie_breaker"]   = NA_REAL,
        _["p.value"]       = pval,
        _["method"]        = "Permutation Mann–Whitney test (B=1)",
        _["alternative"]   = alternative
      );
      if (keep_perm_stats && perm_primary.size() > 0) out["perm_primary"] = perm_primary;
      return out;
    }

    UnpairedKeys K = make_keys_unpaired(x1.size(), y1.size(), B, btype, shift_x, shift_y);
    double S0 = bagged_unpaired_Tdag(x1, y1, raw_statfn, K, B, alternative, ref);

    double pval = perm_bootei_unpaired(x1, y1, raw_statfn,
                                       T0, S0, K, B, R,
                                       alternative, ref,
                                       keep_perm_stats, perm_primary, perm_tie);

    std::string method_str =
      "BOOTEI Mann–Whitney test (tie-breaker via " + describe_boot_type(btype) + ")";

    List out = List::create(
      _["statistic_raw"] = Traw0,
      _["statistic"]     = T0,
      _["tie_breaker"]   = S0,
      _["p.value"]       = pval,
      _["method"]        = method_str,
      _["alternative"]   = alternative
    );
    if (keep_perm_stats) {
      if (perm_primary.size() > 0) out["perm_primary"] = perm_primary;
      if (perm_tie.size() > 0)     out["perm_tie"]     = perm_tie;
    }
    return out;
  }

  // ---------------------------------------------------------------
  // Kruskal–Wallis (paired x + group labels)
  // ---------------------------------------------------------------
  if (test == "kruskalwallis") {
    NumericVector xv(x);
    IntegerVector gy(y);
    if (xv.size() != gy.size())
      stop("For kruskalwallis test x and y must have equal length (response, group).");

    auto raw_statfn = kruskal_factory(gy);

    // For two-sided, ref=0 is natural for H>=0 (still supported for API symmetry)
    const double ref = 0.0;

    double shift = 0.0;
    if (btype == BOOT_SOBOL_SHIFT) {
      RNGScope scope;
      shift = R::runif(0.0, 1.0);
    }

    double Traw0 = raw_statfn(xv, gy);
    double T0    = T_dagger(Traw0, alternative, ref);

    NumericVector perm_primary, perm_tie;

    if (B <= 1) {
      double pval = perm_classical_paired(xv, gy, raw_statfn, T0, R, alternative, ref,
                                          midp, keep_perm_stats, perm_primary);

      List out = List::create(
        _["statistic_raw"] = Traw0,
        _["statistic"]     = T0,
        _["tie_breaker"]   = NA_REAL,
        _["p.value"]       = pval,
        _["method"]        = "Permutation Kruskal–Wallis test (B=1)",
        _["alternative"]   = alternative
      );
      if (keep_perm_stats && perm_primary.size() > 0) out["perm_primary"] = perm_primary;
      return out;
    }

    std::vector<int> keys = make_keys_paired(xv.size(), B, btype, shift);
    double S0 = bagged_paired_Tdag(xv, gy, raw_statfn, keys, B, alternative, ref);

    double pval = perm_bootei_paired(xv, gy, raw_statfn,
                                     T0, S0, keys, B, R,
                                     alternative, ref,
                                     keep_perm_stats, perm_primary, perm_tie);

    std::string method_str =
      "BOOTEI Kruskal–Wallis test (tie-breaker via " + describe_boot_type(btype) + ")";

    List out = List::create(
      _["statistic_raw"] = Traw0,
      _["statistic"]     = T0,
      _["tie_breaker"]   = S0,
      _["p.value"]       = pval,
      _["method"]        = method_str,
      _["alternative"]   = alternative
    );
    if (keep_perm_stats) {
      if (perm_primary.size() > 0) out["perm_primary"] = perm_primary;
      if (perm_tie.size() > 0)     out["perm_tie"]     = perm_tie;
    }
    return out;
  }

  stop("Unknown test type. Use 'chisq', 'spearman', 'mannwhitney', or 'kruskalwallis'.");
  return List::create();
}
