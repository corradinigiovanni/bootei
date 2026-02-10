## ----------------- compile / load -----------------

library(bootei)
library(mlbench)
data("Soybean", package = "mlbench")

dd <- Soybean
dd$Class <- droplevels(dd$Class)          
classes  <- levels(dd$Class)

to_numeric_ordinal <- function(z) {
  if (is.factor(z) || is.character(z)) {
    suppressWarnings(as.numeric(as.character(z)))
  } else {
    suppressWarnings(as.numeric(z))
  }
}

prep_cat <- function(z) {
  if (is.factor(z)) droplevels(z) else droplevels(factor(z))
}

########################################################
######################--- Chi-square ---################
########################################################

sub <- dd[dd$Class == classes[3], , drop = FALSE]

x0 <- sub[[3]]
y0 <- sub[[31]]

ok <- !(is.na(x0) | is.na(y0))
x  <- prep_cat(x0[ok])     
y  <- prep_cat(y0[ok])

tab_xy <- table(x, y)
print(tab_xy)

p_perm_chi2 <- bootei(
  x, y,
  B = 1L, R = 5000L, perm_seed = 910L,
  test = "chisq"
)$p.value

p_bootei_chi2 <- bootei(
  x, y,
  B = 200L, R = 5000L, perm_seed = 910L,
  test = "chisq"
)$p.value

p_perm_chi2 # 0.07258548
p_bootei_chi2 # 0.03379324


######################################################
######################--- t-test ---#############
######################################################


x <- to_numeric_ordinal(dd[dd$Class == classes[3], 9])
y <- to_numeric_ordinal(dd[dd$Class == classes[4], 9])

x <- x[!is.na(x)]
y <- y[!is.na(y)]

c(nx = length(x), ny = length(y))
c(uniq_x = length(unique(x)), uniq_y = length(unique(y)))

p_perm_t <- bootei(
  x, y,
  B = 1L, R = 5000L, perm_seed = 910L,
  test = "welch"
)$p.value

p_bootei_t <- bootei(
  x, y,
  B = 200L, R = 5000L, perm_seed = 910L,
  test = "welch"
)$p.value

p_perm_t # 0.05858828
p_bootei_t # 0.04239152


