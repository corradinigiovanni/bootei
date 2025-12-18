// =====================================================================
//  DETBOOT-EI – Quasi-Random (QMC) Bootstrap Ensemble (Sobol version)
//  (core C++/Rcpp engine, con random shift su Sobol)
//
//  Implementa χ², Spearman, Mann–Whitney e Kruskal–Wallis usando
//  un Efron-bootstrap derandomizzato via sequenze a bassa discrepanza
//  (Sobol 1D, stile Do & Hall 1991 / Aidara 2019), con
//  **random digital shift**:
//
//    u_t      = Sobol1D_raw(t)
//    u'_t     = (u_t + delta) mod 1
//    idx_t    = floor( n * u'_t ), idx_t ∈ {0,...,n-1}.
//
//  Idea QMC:
//
//    - Bootstrap classico = n estrazioni U(0,1) → indici {1,...,n}.
//    - Qui sostituiamo U(0,1) pseudo-casuali con una sequenza
//      quasi-random di Sobol in [0,1), random-shifted:
//
//          u_t_raw  = Sobol1D_raw(t)
//          u_t      = (u_t_raw + shift) mod 1
//
//    - Per il bootstrap b-esimo (b=0,...,B-1) e posizione s=0,...,n-1:
//          t   = b*n + s + 1  (indice globale 1..Bn)
//          u_t = Sobol1D_shifted(t, shift)
//
//    - Unpaired (Mann–Whitney):
//        * per x: uso Sobol1D_shifted(t_x, shift_x)
//        * per y: uso Sobol1D_shifted(t_y + OFFSET, shift_y)
//          (due stream Sobol 1D “quasi-indipendenti” con shift separati).
//
//  Lo statistico ensemble T_bar è la media su T_b = T(resample_b).
//
//  Il p-value è ottenuto via permutation test:
//    per ogni permutazione si ricalcola lo stesso T_bar QMC
//    (stesso random shift; randomizzazione solo sulle etichette).
//
//  Parametri:
//
//    test        = "chisq" | "spearman" | "mannwhitney" | "kruskalwallis"
//    B           = # repliche QMC. Se B <= 1 → stat "grezzo".
//    R           = # permutazioni per il permutation test
//    alternative = "two.sided" | "greater" | "less"
//    perm_seed   = seed per le permutazioni (RNG di R)
//    midp        = se TRUE, usa mid-p-value (solo quando B == 1)
//
//  Compila in R con: Rcpp::sourceCpp("detboot_ei_qmc_sobol.cpp")
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

using namespace Rcpp;

/* ================================================================
 Helper: average ranks (handles ties)
 ================================================================ */
NumericVector rank_numeric(const NumericVector &x) {
  int n = x.size();
  std::vector< std::pair<double,int> > ord(n);
  for (int i = 0; i < n; ++i) ord[i] = std::make_pair(x[i], i);
  std::sort(ord.begin(), ord.end(),
            [](const std::pair<double,int>& a,
               const std::pair<double,int>& b){ return a.first < b.first; });
  NumericVector rk(n);
  int i = 0;
  while (i < n) {
    int j = i;
    while (j + 1 < n && ord[j+1].first == ord[i].first) ++j;
    double r = (i + j + 2) / 2.0;  // 1-based average rank
    for (int k = i; k <= j; ++k) rk[ ord[k].second ] = r;
    i = j + 1;
  }
  return rk;
}

/* small helper per NA double */
inline bool is_na_double(double x) {
  return NumericVector::is_na(x);
}

/* ================================================================
 χ² for independence – versione ottimizzata:
 - comprimiamo x,y in indici di riga/colonna interi
 - bootstrap QMC e permutazioni lavorano su (row_idx, col_idx)
 senza toccare CharacterVector a ogni resample
 ================================================================ */

inline double chi2_from_table(const std::vector< std::vector<double> > &tab,
                              const std::vector<double> &rs,
                              const std::vector<double> &cs,
                              double tot){
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

// χ² per il campione osservato a partire da (row_idx, col_idx)
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

/* ================================================================
 Kruskal–Wallis test (multi-sample, ordinal/continuous)
 - x: NumericVector (response)
 - g: IntegerVector (group labels)
 - factory fissa i livelli di g per evitare problemi nei resample
 ================================================================ */

std::function<double(const NumericVector&, const IntegerVector&)>
kruskal_factory(const IntegerVector &g_ref) {
  std::unordered_set<int> seen;
  std::vector<int> levels;
  levels.reserve(g_ref.size());
  
  for (int i = 0; i < g_ref.size(); ++i) {
    if (IntegerVector::is_na(g_ref[i])) continue;
    int v = g_ref[i];
    if (!seen.count(v)) {
      seen.insert(v);
      levels.push_back(v);
    }
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
    
    // Compatta rimuovendo NA e livelli sconosciuti
    std::vector<double> x2;
    std::vector<int> grp2;
    x2.reserve(nloc);
    grp2.reserve(nloc);
    
    for (int i = 0; i < nloc; ++i) {
      if (NumericVector::is_na(x[i]) || IntegerVector::is_na(g[i])) continue;
      auto it = gmap.find(g[i]);
      if (it == gmap.end()) continue;
      x2.push_back(x[i]);
      grp2.push_back(it->second); // 0..K-1
    }
    
    int N = (int)x2.size();
    if (N < 2) return NA_REAL;
    
    // Ranks sui dati poolati
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
    for (int j = 0; j < K; ++j)
      if (nj[j] > 0) k_eff++;
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

/* ================================================================
 Spearman ρ and Mann–Whitney U
 ================================================================ */

double spearman_stat(const NumericVector &x, const NumericVector &y){
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

double mannwhitney_U(const NumericVector &x, const NumericVector &y){
  int nx = x.size(), ny = y.size();
  if (nx == 0 || ny == 0) return NA_REAL;
  
  int n = nx + ny;
  NumericVector z(n);
  IntegerVector grp(n);
  for (int i = 0; i < nx; ++i) { z[i] = x[i];       grp[i]     = 0; }
  for (int j = 0; j < ny; ++j) { z[nx + j] = y[j];  grp[nx+j] = 1; }
  
  NumericVector rk = rank_numeric(z);
  double R1 = 0.0;
  for (int i = 0; i < n; ++i)
    if (grp[i] == 0) R1 += rk[i];
    
    return R1 - (double)nx * (nx + 1) / 2.0;
}

/* ================================================================
 Sobol 1D (base 2) – direct i-th point via Gray code + random shift
 ================================================================ */

// Sobol grezzo in [0,1) per indice 1-based 'index'
inline double sobol_1d_raw(uint32_t index) {
  if (index == 0u) return 0.0;
  const int BITS = 32;
  uint32_t X = 0u;
  uint32_t g = index ^ (index >> 1); // Gray code
  // direction numbers: v_j = 1 / 2^j, rappresentati come int
  for (int k = 0; k < BITS; ++k) {
    if (g & (1u << k)) {
      // XOR con v_{k+1} = 1 << (BITS-1-k)
      X ^= (1u << (BITS - 1 - k));
    }
  }
  const double norm = 4294967296.0; // 2^32
  return (double)X / norm;
}

// versione con random digital shift in [0,1)
inline double sobol_1d(uint32_t index, double shift) {
  double u = sobol_1d_raw(index);
  u += shift;
  if (u >= 1.0) u -= 1.0;
  return u;
}

/* ================================================================
 QMC bootstrap helpers (Sobol + random shift)
 ================================================================ */

/*
 Schema paired (x[i], y[i] trattati come coppie):
 
 - n osservazioni, B repliche.
 - Se B <= 1: stat grezzo.
 - Per b=0..B-1, s=0..n-1:
 t       = b*n + s + 1   (1-based)
 u'_t    = Sobol1D(t, shift)
 idx     = floor( n * u'_t )
 xb[s]   = x[idx]; yb[s] = y[idx]
 */

template<typename VECX, typename VECY, typename STATFN>
double qmc_paired(const VECX &x, const VECY &y,
                  STATFN statfn, int B, double shift) {
  int n = x.size();
  if (n != (int)y.size()) return NA_REAL;
  if (n <= 1 || B <= 1) return statfn(x, y);
  
  VECX xb(n);
  VECY yb(n);
  
  double acc = 0.0;
  int ok = 0;
  
  for (int b = 0; b < B; ++b) {
    for (int s = 0; s < n; ++s) {
      long long t_ll = (long long)b * (long long)n + (long long)s + 1LL;
      uint32_t t = (uint32_t)(t_ll & 0xFFFFFFFFULL); // wrap mod 2^32
      double u = sobol_1d(t, shift);
      int idx = (int)std::floor((double)n * u);
      if (idx < 0)   idx = 0;
      if (idx >= n)  idx = n - 1;
      xb[s] = x[idx];
      yb[s] = y[idx];
    }
    double val = statfn(xb, yb);
    if (!is_na_double(val)) {
      acc += val;
      ++ok;
    }
  }
  
  if (ok == 0) return statfn(x, y);
  return acc / (double)ok;
}

/*
 Schema unpaired (Mann–Whitney):
 
 - x di lunghezza nx, y di lunghezza ny.
 - Se B <= 1: stat grezzo.
 - Per ciascun b:
 per x: u_x = Sobol1D(t_x, shift_x)
 per y: u_y = Sobol1D(t_y + OFFSET, shift_y)
 */

double qmc_unpaired(
    const NumericVector &x,
    const NumericVector &y,
    const std::function<double(const NumericVector&, const NumericVector&)> &statfn,
    int B,
    double shift_x,
    double shift_y) {
  
  int nx = x.size(), ny = y.size();
  if (nx <= 0 || ny <= 0 || B <= 1) return statfn(x, y);
  
  NumericVector xb(nx), yb(ny);
  
  double acc = 0.0;
  int ok = 0;
  const uint32_t OFFSET = 123456789u; // offset per secondo stream Sobol
  
  for (int b = 0; b < B; ++b) {
    // gruppo x con Sobol1D(t, shift_x)
    for (int s = 0; s < nx; ++s) {
      long long t_ll = (long long)b * (long long)nx + (long long)s + 1LL;
      uint32_t t = (uint32_t)(t_ll & 0xFFFFFFFFULL);
      double u = sobol_1d(t, shift_x);
      int idx = (int)std::floor((double)nx * u);
      if (idx < 0)    idx = 0;
      if (idx >= nx)  idx = nx - 1;
      xb[s] = x[idx];
    }
    // gruppo y con Sobol1D(t + OFFSET, shift_y)
    for (int s = 0; s < ny; ++s) {
      long long t_ll = (long long)b * (long long)ny + (long long)s + 1LL;
      uint32_t t = (uint32_t)(t_ll & 0xFFFFFFFFULL);
      uint32_t t2 = t + OFFSET;
      double u = sobol_1d(t2, shift_y);
      int idx = (int)std::floor((double)ny * u);
      if (idx < 0)    idx = 0;
      if (idx >= ny)  idx = ny - 1;
      yb[s] = y[idx];
    }
    
    double val = statfn(xb, yb);
    if (!is_na_double(val)) {
      acc += val;
      ++ok;
    }
  }
  
  if (ok == 0) return statfn(x, y);
  return acc / (double)ok;
}

/* ================================================================
 QMC per χ² con (row_idx, col_idx) interi
 ================================================================ */

// QMC bootstrap paired per χ² (lavorando sui soli indici di riga/colonna)
double qmc_chisq_int(const IntegerVector &row_idx,
                     const IntegerVector &col_idx,
                     int nr, int nc,
                     int B,
                     double shift) {
  int n = row_idx.size();
  if (n != col_idx.size()) return NA_REAL;
  if (n <= 1 || B <= 1) return chisq_from_pairs_int(row_idx, col_idx, nr, nc);
  
  double acc = 0.0;
  int ok = 0;
  
  // vettori di contatori allocati una volta sola e azzerati per ogni b
  std::vector< std::vector<double> > tab(nr, std::vector<double>(nc, 0.0));
  std::vector<double> rs(nr, 0.0), cs(nc, 0.0);
  
  for (int b = 0; b < B; ++b) {
    // azzera tabella e margini
    for (int r = 0; r < nr; ++r) {
      std::fill(tab[r].begin(), tab[r].end(), 0.0);
      rs[r] = 0.0;
    }
    std::fill(cs.begin(), cs.end(), 0.0);
    double tot = 0.0;
    
    for (int s = 0; s < n; ++s) {
      long long t_ll = (long long)b * (long long)n + (long long)s + 1LL;
      uint32_t t = (uint32_t)(t_ll & 0xFFFFFFFFULL);
      double u = sobol_1d(t, shift);
      int idx = (int)std::floor((double)n * u);
      if (idx < 0)   idx = 0;
      if (idx >= n)  idx = n - 1;
      
      int r = row_idx[idx];
      int c = col_idx[idx];
      if (r < 0 || r >= nr || c < 0 || c >= nc) continue;
      
      tab[r][c] += 1.0;
      rs[r]     += 1.0;
      cs[c]     += 1.0;
      tot       += 1.0;
    }
    
    double val = chi2_from_table(tab, rs, cs, tot);
    if (!is_na_double(val)) {
      acc += val;
      ++ok;
    }
  }
  
  if (ok == 0) return chisq_from_pairs_int(row_idx, col_idx, nr, nc);
  return acc / (double)ok;
}

/* ================================================================
 Permutation helpers
 - usano qmc_* come statistica ensemble (con shift fisso)
 - seed delle permutazioni passato da bootei tramite set.seed
 - midp: se TRUE e B == 1 → usa mid-p senza add-one
 ================================================================ */

template<typename VECX, typename VECY, typename STATFN>
double perm_paired(const VECX &x, const VECY &y, STATFN statfn,
                   int B, int R, const std::string &alt,
                   double ref, bool midp,
                   double shift) {
  
  double obs = qmc_paired(x, y, statfn, B, shift);
  if (R <= 0 || is_na_double(obs)) return NA_REAL;
  
  RNGScope scope;
  int n = x.size();
  
  int ge_or_eq = 0;
  int gt = 0, eq = 0;
  int valid = 0;
  const double eps = 1e-12;
  
  for (int r = 0; r < R; ++r) {
    IntegerVector perm = Rcpp::sample(n, n, false);
    VECY yp(n);
    for (int i = 0; i < n; ++i)
      yp[i] = y[perm[i] - 1];
    
    double s = qmc_paired(x, yp, statfn, B, shift);
    if (is_na_double(s)) continue;
    ++valid;
    
    if (alt == "greater") {
      if (s >= obs - eps) ge_or_eq++;
      if (s >  obs + eps) gt++;
      else if (std::fabs(s - obs) <= eps) eq++;
    } else if (alt == "less") {
      if (s <= obs + eps) ge_or_eq++;
      if (s <  obs - eps) gt++;
      else if (std::fabs(s - obs) <= eps) eq++;
    } else { // two.sided
      double dobj = std::fabs(obs - ref);
      double d = std::fabs(s - ref);
      if (d >= dobj - eps) ge_or_eq++;
      if (d >  dobj + eps) gt++;
      else if (std::fabs(d - dobj) <= eps) eq++;
    }
  }
  
  if (valid == 0) return NA_REAL;
  
  if (midp && B == 1) {
    return ((double)gt + 0.5 * (double)eq) / (double)valid;
  } else {
    return (double)(ge_or_eq + 1) / (double)(valid + 1);
  }
}

double perm_unpaired(const NumericVector &x, const NumericVector &y,
                     const std::function<double(const NumericVector&, const NumericVector&)> &statfn,
                     int B, int R, const std::string &alt,
                     double ref, bool midp,
                     double shift_x,
                     double shift_y) {
  
  double obs = qmc_unpaired(x, y, statfn, B, shift_x, shift_y);
  if (R <= 0 || is_na_double(obs)) return NA_REAL;
  
  RNGScope scope;
  int nx = x.size(), ny = y.size();
  int N  = nx + ny;
  if (N == 0) return NA_REAL;
  
  NumericVector pool(N);
  for (int i = 0; i < nx; ++i) pool[i]    = x[i];
  for (int j = 0; j < ny; ++j) pool[nx+j] = y[j];
  
  int ge_or_eq = 0;
  int gt = 0, eq = 0;
  int valid = 0;
  const double eps = 1e-12;
  
  for (int r = 0; r < R; ++r) {
    IntegerVector perm = Rcpp::sample(N, N, false);
    NumericVector xg(nx), yg(ny);
    for (int i = 0; i < nx; ++i) xg[i] = pool[perm[i]      - 1];
    for (int j = 0; j < ny; ++j) yg[j] = pool[perm[nx + j] - 1];
    
    double s = qmc_unpaired(xg, yg, statfn, B, shift_x, shift_y);
    if (is_na_double(s)) continue;
    ++valid;
    
    if (alt == "greater") {
      if (s >= obs - eps) ge_or_eq++;
      if (s >  obs + eps) gt++;
      else if (std::fabs(s - obs) <= eps) eq++;
    } else if (alt == "less") {
      if (s <= obs + eps) ge_or_eq++;
      if (s <  obs - eps) gt++;
      else if (std::fabs(s - obs) <= eps) eq++;
    } else { // two.sided
      double dobj = std::fabs(obs - ref);
      double d = std::fabs(s - ref);
      if (d >= dobj - eps) ge_or_eq++;
      if (d >  dobj + eps) gt++;
      else if (std::fabs(d - dobj) <= eps) eq++;
    }
  }
  
  if (valid == 0) return NA_REAL;
  
  if (midp && B == 1) {
    return ((double)gt + 0.5 * (double)eq) / (double)valid;
  } else {
    return (double)(ge_or_eq + 1) / (double)(valid + 1);
  }
}

/* ================================================================
 Permutation helper specializzato per χ² con row/col int
 ================================================================ */

double perm_chisq_int(const IntegerVector &row_idx,
                      const IntegerVector &col_idx,
                      int nr, int nc,
                      int B, int R,
                      bool midp,
                      double shift) {
  
  double obs = qmc_chisq_int(row_idx, col_idx, nr, nc, B, shift);
  if (R <= 0 || is_na_double(obs)) return NA_REAL;
  
  RNGScope scope;
  int n = row_idx.size();
  if (n != col_idx.size()) return NA_REAL;
  
  int ge_or_eq = 0;
  int gt = 0, eq = 0;
  int valid = 0;
  const double eps = 1e-12;
  
  // vettori temporanei per col_idx permutato
  IntegerVector col_perm(n);
  
  for (int r = 0; r < R; ++r) {
    IntegerVector perm = Rcpp::sample(n, n, false);
    for (int i = 0; i < n; ++i) {
      col_perm[i] = col_idx[ perm[i] - 1 ];
    }
    
    double s = qmc_chisq_int(row_idx, col_perm, nr, nc, B, shift);
    if (is_na_double(s)) continue;
    ++valid;
    
    // χ² è sempre greater-tailed
    if (s >= obs - eps) ge_or_eq++;
    if (s >  obs + eps) gt++;
    else if (std::fabs(s - obs) <= eps) eq++;
  }
  
  if (valid == 0) return NA_REAL;
  
  if (midp && B == 1) {
    return ((double)gt + 0.5 * (double)eq) / (double)valid;
  } else {
    return (double)(ge_or_eq + 1) / (double)(valid + 1);
  }
}

/* ================================================================
 Exported wrapper
 ================================================================ */

// [[Rcpp::export]]
List bootei(SEXP x, SEXP y,
            std::string test = "chisq",
            int B = 100,
            int R = 1000,
            std::string alternative = "two.sided",
            double perm_seed = NA_REAL,
            bool midp = false) {
  
  if (alternative != "two.sided" &&
      alternative != "greater" &&
      alternative != "less") {
    stop("alternative must be 'two.sided', 'greater', or 'less'");
  }
  
  // Imposta il seed per le permutazioni (se fornito)
  if (!NumericVector::is_na(perm_seed)) {
    Environment base = Environment::namespace_env("base");
    Function set_seed = base["set.seed"];
    set_seed(perm_seed);
  }
  
  if (test == "chisq") {
    CharacterVector xc(x), yc(y);
    if (xc.size() != yc.size())
      stop("For chisq test x and y must have equal length.");
    
    int n = xc.size();
    if (n <= 1) {
      return List::create(
        _["statistic"]    = NA_REAL,
        _["p.value"]      = NA_REAL,
        _["method"]       = "DETBOOT-EI χ² test (Sobol QMC bootstrap ensemble)",
        _["alternative"]  = "greater"
      );
    }
    
    // mappiamo i livelli di riga/colonna in interi 0..nr-1, 0..nc-1
    std::unordered_map<std::string,int> rmap, cmap;
    rmap.reserve(n);
    cmap.reserve(n);
    
    for (int i = 0; i < n; ++i) {
      if (CharacterVector::is_na(xc[i]) || CharacterVector::is_na(yc[i])) continue;
      rmap[ as<std::string>(xc[i]) ];
      cmap[ as<std::string>(yc[i]) ];
    }
    
    int id = 0;
    for (auto &kv : rmap) kv.second = id++;
    int nr = id;
    
    id = 0;
    for (auto &kv : cmap) kv.second = id++;
    int nc = id;
    
    if (nr < 2 || nc < 2) {
      return List::create(
        _["statistic"]    = NA_REAL,
        _["p.value"]      = NA_REAL,
        _["method"]       = "DETBOOT-EI χ² test (Sobol QMC bootstrap ensemble)",
        _["alternative"]  = "greater"
      );
    }
    
    IntegerVector row_idx(n), col_idx(n);
    for (int i = 0; i < n; ++i) {
      if (CharacterVector::is_na(xc[i]) || CharacterVector::is_na(yc[i])) {
        row_idx[i] = -1;
        col_idx[i] = -1;
      } else {
        auto it_r = rmap.find( as<std::string>(xc[i]) );
        auto it_c = cmap.find( as<std::string>(yc[i]) );
        if (it_r == rmap.end() || it_c == cmap.end()) {
          row_idx[i] = -1;
          col_idx[i] = -1;
        } else {
          row_idx[i] = it_r->second;
          col_idx[i] = it_c->second;
        }
      }
    }
    
    // random shift QMC (digital shift) deterministico dato perm_seed
    {
      RNGScope scope;
      double shift = R::runif(0.0, 1.0);
      
      double stat = qmc_chisq_int(row_idx, col_idx, nr, nc, B, shift);
      double pval = perm_chisq_int(row_idx, col_idx, nr, nc, B, R, midp, shift);
      
      return List::create(
        _["statistic"]    = stat,
        _["p.value"]      = pval,
        _["method"]       = "DETBOOT-EI χ² test (Sobol QMC bootstrap ensemble, random-shift Sobol) ",
        _["alternative"]  = "greater"
      );
    }
  }
  
  else if (test == "spearman") {
    NumericVector xn(x), yn(y);
    if (xn.size() != yn.size())
      stop("For spearman test x and y must have equal length.");
    
    auto statfn = [](const NumericVector &a, const NumericVector &b){
      return spearman_stat(a, b);
    };
    
    RNGScope scope;
    double shift = R::runif(0.0, 1.0);
    
    double stat = qmc_paired(xn, yn, statfn, B, shift);
    double pval = perm_paired(xn, yn, statfn, B, R, alternative, 0.0, midp, shift);
    
    return List::create(
      _["statistic"]    = stat,
      _["p.value"]      = pval,
      _["method"]       = "DETBOOT-EI Spearman test (Sobol QMC bootstrap ensemble, random-shift)",
      _["alternative"]  = alternative
    );
  }
  
  else if (test == "mannwhitney") {
    NumericVector x1(x), y1(y);
    
    auto statfn = [](const NumericVector &a, const NumericVector &b){
      return mannwhitney_U(a, b);
    };
    
    RNGScope scope;
    double shift_x = R::runif(0.0, 1.0);
    double shift_y = R::runif(0.0, 1.0);
    
    double stat = qmc_unpaired(x1, y1, statfn, B, shift_x, shift_y);
    double ref  = (double)x1.size() * (double)y1.size() / 2.0;
    double pval = perm_unpaired(x1, y1, statfn, B, R, alternative, ref, midp,
                                shift_x, shift_y);
    
    return List::create(
      _["statistic"]    = stat,
      _["p.value"]      = pval,
      _["method"]       = "DETBOOT-EI Mann–Whitney test (Sobol QMC bootstrap ensemble, random-shift)",
      _["alternative"]  = alternative
    );
  }
  
  else if (test == "kruskalwallis") {
    NumericVector xv(x);
    IntegerVector gy(y);
    if (xv.size() != gy.size())
      stop("For kruskalwallis test x and y must have equal length (response, group).");
    
    auto statfn = kruskal_factory(gy);
    
    RNGScope scope;
    double shift = R::runif(0.0, 1.0);
    
    double stat = qmc_paired(xv, gy, statfn, B, shift);
    double ref  = 0.0; // H ≥ 0, ref per two.sided
    double pval = perm_paired(xv, gy, statfn, B, R, alternative, ref, midp, shift);
    
    return List::create(
      _["statistic"]    = stat,
      _["p.value"]      = pval,
      _["method"]       = "DETBOOT-EI Kruskal–Wallis test (Sobol QMC bootstrap ensemble, random-shift)",
      _["alternative"]  = alternative
    );
  }
  
  else {
    stop("Unknown test type. Use 'chisq', 'spearman', 'mannwhitney', or 'kruskalwallis'.");
  }
}
