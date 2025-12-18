
Rcpp::sourceCpp('C:/Users/innav/Desktop/bootei/bootei.cpp')

library(mlbench)
data("Soybean")

dd <- Soybean
classes <- unique(Soybean$Class)

# Spearman

res <- list()
k <- 0

for(cl1 in 1:19){
  for(i in c(2,4,5,9,11)){
    for(j in c(2,4,5,9,11)){
      if((i != j)){
        seme <- i + j + cl1 
        bb <- bootei(dd[dd$Class == classes[cl1],i], dd[dd$Class == classes[cl1],j], B = 1,
                     perm_seed = seme, test = 'spearman')$p.value
        if(!is.na(bb)){
          if(bb > 0.05 & bb < 0.2){
            
            bb2 <- bootei(dd[dd$Class == classes[cl1],i], dd[dd$Class == classes[cl1],j], B = 100,
                          perm_seed = seme, test = 'spearman')$p.value
            if(bb2 < 0.05){
              k <- k + 1
              res[[k]] <- c(cl1, i, j, bb, bb2)
              
            }
          }
        }
      }
    }
  }
}


resu <- matrix(0, nrow = 2, ncol = 10)

for(i in 1:10){
  
  print(i)
  
  seme <- 910 - i
  
  resu[1,i] <- (bootei(dd[dd$Class == classes[3],9], 
                       dd[dd$Class == classes[3],11], 
                       B = 1, R = 1e4, 
                       perm_seed = seme, test = 'spearman')$p.value > 0.05) &
    (bootei(dd[dd$Class == classes[3],9], 
            dd[dd$Class == classes[3],11], 
            B = 5e2, R = 1e4, 
            perm_seed = seme, test = 'spearman')$p.value < 0.05)
  
  
  resu[2,i] <- (bootei(dd[dd$Class == classes[8],4], 
                       dd[dd$Class == classes[8],11], 
                       B = 1, R = 1e4, 
                       perm_seed = seme, test = 'spearman')$p.value > 0.05) &
    (bootei(dd[dd$Class == classes[8],4], 
            dd[dd$Class == classes[8],11], 
            B = 5e2, R = 1e4, 
            perm_seed = seme, test = 'spearman')$p.value < 0.05)
  
  
  
}

rowSums(resu)/10 # 1 1

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

# Mann Whitney

res <- list()
k <- 0

for(cl1 in 1:19){
  for(cl2 in 1:19){
   for(i in c(2,4,5,9,11)){
     for(j in c(2,4,5,9,11)){
       
        seme <- i + j + cl1 
        bb <- bootei(dd[dd$Class == classes[cl1],i], dd[dd$Class == classes[cl2],j], B = 1, R = 1e2,
                     perm_seed = seme, test = 'mannwhitney')$p.value
        if(!is.na(bb)){
          if(bb > 0.05 & bb < 0.2){
            
            bb2 <- bootei(dd[dd$Class == classes[cl1],i], dd[dd$Class == classes[cl2],j], B = 100, R = 1e2,
                          perm_seed = seme, test = 'mannwhitney')$p.value
            if(bb2 < 0.05){
              k <- k + 1
              res[[k]] <- c(cl1, cl2, i, j, bb, bb2)
              
            }
          }
        }
      }
    }
  }
}





resu <- matrix(0, nrow = 4, ncol = 10)

for(i in 1:10){
  
  print(i)
  
  seme <- 910 * i
  
  
  
  resu[1,i] <- (bootei(dd[dd$Class == classes[4],5], 
                       dd[dd$Class == classes[17],4], 
                       B = 1, R = 1e4, 
                       perm_seed = seme, test = 'mannwhitney')$p.value > 0.05) &
    (bootei(dd[dd$Class == classes[4],5], 
            dd[dd$Class == classes[17],4], 
            B = 5e2, R = 1e4, 
            perm_seed = seme, test = 'mannwhitney')$p.value < 0.05)
  
  
  resu[2,i] <- (bootei(dd[dd$Class == classes[5],2], 
                       dd[dd$Class == classes[17],4], 
                       B = 1, R = 1e4, 
                       perm_seed = seme, test = 'mannwhitney')$p.value > 0.05) &
    (bootei(dd[dd$Class == classes[5],2], 
            dd[dd$Class == classes[17],4], 
            B = 5e2, R = 1e4, 
            perm_seed = seme, test = 'mannwhitney')$p.value < 0.05)
  
  resu[3,i] <- (bootei(dd[dd$Class == classes[5],2], 
                       dd[dd$Class == classes[17],5], 
                       B = 1, R = 1e4, 
                       perm_seed = seme, test = 'mannwhitney')$p.value > 0.05) &
    (bootei(dd[dd$Class == classes[5],2], 
            dd[dd$Class == classes[17],5], 
            B = 5e2, R = 1e4, 
            perm_seed = seme, test = 'mannwhitney')$p.value < 0.05)
  
  resu[4,i] <- (bootei(dd[dd$Class == classes[5],11], 
                       dd[dd$Class == classes[17],5], 
                       B = 1, R = 1e4, 
                       perm_seed = seme, test = 'mannwhitney')$p.value > 0.05) &
    (bootei(dd[dd$Class == classes[5],11], 
            dd[dd$Class == classes[17],5], 
            B = 5e2, R = 1e4, 
            perm_seed = seme, test = 'mannwhitney')$p.value < 0.05)
  
  
  
}

rowSums(resu)/10


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


# Chi2

res <- list()
k <- 0

for(cl1 in 1:19){
  for(i in c(3,6:36)){
    for(j in c(3,6:36)){
      if((i != j)){
        seme <- i + j + cl1 
        bb <- bootei(dd[dd$Class == classes[cl1],i], dd[dd$Class == classes[cl1],j], B = 1,
                     perm_seed = seme, test = 'chisq')$p.value
        if(!is.na(bb)){
          if(bb > 0.05 & bb < 0.4){
            
            bb2 <- bootei(dd[dd$Class == classes[cl1],i], dd[dd$Class == classes[cl1],j], B = 100,
                          perm_seed = seme, test = 'chisq')$p.value
            if(bb2 < 0.05){
              k <- k + 1
              res[[k]] <- c(cl1, i, j, bb, bb2)
              
            }
          }
        }
      }
    }
  }
}


resu <- numeric(10)

for(i in 1:10){
  
  print(i)
  
  seme <- i
  
  resu[1,i] <- (bootei(dd[dd$Class == classes[15],29], 
                       dd[dd$Class == classes[15],12], 
                       B = 1, R = 1e4, 
                       perm_seed = seme, test = 'chisq')$p.value > 0.05) &
    (bootei(dd[dd$Class == classes[15],29], 
            dd[dd$Class == classes[15],12], 
            B = 1e3, R = 1e4, 
            perm_seed = seme, test = 'chisq')$p.value < 0.05)
  
}

sum(resu)/10 # 1


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

## ===============================================================
## Soybean search loop: Kruskal–Wallis with different classes (cl1, cl2, cl3)
## Goal: find cases where BOOTEI rejects (p < 0.05) but classical permutation KW does not.
##
## Construction (3-sample KW):
##   - pick three (not necessarily distinct) classes cl1, cl2, cl3
##   - pick ONE attribute column i used as ordinal response for all 3 groups
##   - define groups by class membership (3 groups)
##   - run KW on pooled data (x = response, g = group labels 1/2/3)
##
## Speed filters:
##   - drop NA
##   - each group must have >= n_min_per_group observations after NA-drop
##   - response must have >= 3 unique values overall (avoid binary-only)
##   - within pooled sample, at least one group must have >= 3 unique values
##     (so we do not only test 3 binary samples)
##   - optional pre-screen: p_class in (0.05, p_upper) before running BOOTEI
##
## Requires: bootei.cpp compiled via Rcpp::sourceCpp
## ===============================================================

## ===============================================================
## Soybean search loop: Kruskal–Wallis with 3, 4, or 5 groups (classes)
## Goal: find cases where BOOTEI rejects (p < 0.05) but classical permutation KW does not.
##
## Groups are defined by selecting K distinct Class levels (K=3,4,5),
## pooling the chosen classes, and running KW on one attribute column.
##
## Speed filters:
##   - drop NA within each group
##   - each group must have >= n_min_per_group obs
##   - response has >= 3 unique overall
##   - at least one group has >= 3 unique values
##   - optional pre-screen: p_class in (0.05, p_upper) before running BOOTEI
## ===============================================================

Rcpp::sourceCpp("C:/Users/innav/Desktop/bootstrap ensemble inference/bootei.cpp")

library(mlbench)
data("Soybean", package = "mlbench")

dd <- Soybean
classes <- unique(dd$Class)
G <- length(classes)

## ---------------- helpers ----------------
is_constant <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) == 0) return(TRUE)
  length(unique(v)) <= 1
}

## Build pooled (x,g) from a vector of class indices and one attribute column
build_kw_data <- function(class_idx_vec, attr_col, dd, classes) {
  xs <- vector("list", length(class_idx_vec))
  ns <- integer(length(class_idx_vec))
  uxs <- integer(length(class_idx_vec))
  
  for (k in seq_along(class_idx_vec)) {
    cl <- class_idx_vec[k]
    xk <- dd[dd$Class == classes[cl], attr_col]
    xk <- xk[!is.na(xk)]
    xs[[k]] <- xk
    ns[k] <- length(xk)
    uxs[k] <- length(unique(xk))
  }
  
  x <- do.call(c, xs)
  g <- unlist(mapply(function(k, n) rep.int(k, n),
                     k = seq_along(ns), n = ns, SIMPLIFY = FALSE),
              use.names = FALSE)
  
  list(x = x, g = as.integer(g), ns = ns, uxs = uxs,
       uniq_overall = length(unique(x)))
}

## Deterministic seed from class set + attribute
seed_from <- function(class_idx_vec, attr_col) {
  v <- sort(as.integer(class_idx_vec))
  s <- (sum(v * c(1L, 97L, 997L, 9973L, 99991L)[seq_along(v)]) + 1000003L * as.integer(attr_col)) %% 2147483647L
  s <- as.integer(s)
  if (is.na(s) || s <= 0L) s <- 1L
  s
}


## ---------------- parameters ----------------
alpha   <- 0.05
R_perm  <- 200L     # screening R (increase later)
B_big   <- 200L
p_upper <- 0.20

attr_cols <- c(2,4,5,9,11)   # expand if desired: 2:36
n_min_per_group <- 5L

## K values to scan
K_set <- 3:10

res <- list()
k_hits <- 0L

## ---------------------------------------------------------------
## Loop over K=3,4,5 and all combinations of distinct classes
## ---------------------------------------------------------------
for (K in K_set) {
  
  combs <- combn(G, K)  # each column is a combination of class indices
  
  for (cc in seq_len(ncol(combs))) {
    
    class_idx_vec <- combs[, cc]
    
    ## Quick size screen before attributes
    n_by_class <- vapply(class_idx_vec,
                         function(cl) sum(dd$Class == classes[cl]),
                         numeric(1))
    if (any(n_by_class < n_min_per_group)) next
    
    for (i in attr_cols) {
      
      dat <- build_kw_data(class_idx_vec, i, dd, classes)
      
      ## Group-size filter after NA-drop
      if (any(dat$ns < n_min_per_group)) next
      
      ## Avoid totally constant across all groups
      if (dat$uniq_overall < 2) next
      
      ## Response must have >= 3 unique overall
      if (dat$uniq_overall < 3) next
      
      ## At least one group has >= 3 unique values
      if (max(dat$uxs) < 3) next
      
      ## Deterministic seed
      seme <- seed_from(class_idx_vec, i)
      
      ## Classical KW permutation p-value (B=1)
      p_class <- bootei(
        dat$x, dat$g,
        test        = "kruskalwallis",
        B           = 1L,
        R           = R_perm,
        alternative = "greater",
        perm_seed   = seme,
        boot_type   = "sobol"
      )$p.value
      
      if (is.na(p_class)) next
      if (!(p_class > alpha && p_class < p_upper)) next
      
      ## BOOTEI KW p-value
      p_bt <- bootei(
        dat$x, dat$g,
        test        = "kruskalwallis",
        B           = B_big,
        R           = R_perm,
        alternative = "greater",
        perm_seed   = seme,
        boot_type   = "sobol"
      )$p.value
      
      if (is.na(p_bt)) next
      
      if (p_bt < alpha) {
        k_hits <- k_hits + 1L
        
        res[[k_hits]] <- data.frame(
          K = K,
          class_idx = paste(class_idx_vec, collapse = ","),
          class_names = paste(as.character(classes[class_idx_vec]), collapse = " | "),
          attr_col  = i,
          attr_name = names(dd)[i],
          ns = paste(dat$ns, collapse = ","),
          uniqs = paste(dat$uxs, collapse = ","),
          uniq_overall = dat$uniq_overall,
          p_classical = p_class,
          p_bootei    = p_bt,
          seed = seme,
          stringsAsFactors = FALSE
        )
        
        cat(sprintf("HIT %d: K=%d classes={%s} attr=%s | p_class=%.4f -> p_bt=%.4f\n",
                    k_hits, K, paste(class_idx_vec, collapse = ","),
                    names(dd)[i], p_class, p_bt))
      }
      
    } # attr
  } # combination
} # K

res_df <- if (length(res) > 0) do.call(rbind, res) else data.frame()
res_df <- res_df[order(res_df$K, res_df$p_bootei, res_df$p_classical), , drop = FALSE]
res_df

## write.csv(res_df, "kw_bootei_hits_K3to5_screening.csv", row.names = FALSE)

## ---------------------------------------------------------------
## Optional: confirmation step for one hit (larger R/B)
## ---------------------------------------------------------------
confirm_kw_hit <- function(hit_row, R_confirm = 10000L, B_confirm = 500L,
                           boot_type = "sobol") {
  
  class_idx_vec <- as.integer(strsplit(hit_row$class_idx, ",")[[1]])
  i <- hit_row$attr_col
  seme <- hit_row$seed
  
  dat <- build_kw_data(class_idx_vec, i, dd, classes)
  
  p_class <- bootei(
    dat$x, dat$g,
    test="kruskalwallis",
    B=1L, R=R_confirm, alternative="greater",
    perm_seed=seme, boot_type=boot_type
  )$p.value
  
  p_bt <- bootei(
    dat$x, dat$g,
    test="kruskalwallis",
    B=B_confirm, R=R_confirm, alternative="greater",
    perm_seed=seme, boot_type=boot_type
  )$p.value
  
  c(p_classical_confirm = p_class, p_bootei_confirm = p_bt)
}

## Example:
## if (nrow(res_df) > 0) confirm_kw_hit(res_df[1, ], R_confirm=10000L, B_confirm=500L)
