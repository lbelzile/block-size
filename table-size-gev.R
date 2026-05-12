library(simsalapar)
setwd("/home/lbelzile/Downloads/blocksize")
out_weibull <- readRDS("power-study-blocksize-1_2.rds")
out_normal <- readRDS("power-study-blocksize-1_3.rds")
B <- 2000L
delta_seq <- seq(1, 18, by = 0.25)
nobs_seq <- c(25L, 50L, 100L)
m_seq <- c(2L, 5L, 10L)

# Run simulations
#Set list of variables for the simulation study
varList <- simsalapar::varlist(
  n.sim = list(type = "N", expr = quote(N), value = B),
  delta = list(type = "grid", value = delta_seq),
  nobs = list(type = "grid", value = nobs_seq),
  m = list(type = "grid", value = m_seq),
  alt = list(type = "inner", value = 1:3)
)
mk <- simsalapar::mkAL(x = out_weibull, vList = varList, repFirst = TRUE)
resArray <- getArray(mk, comp = "value")

# Compute percentage of rejection for a 5% test under H0
# Kolmogorov-Smirnov tests of uniformity
size_ar_weibull <- apply(
  X = resArray[, 1, , , ],
  MARGIN = 1:3,
  FUN = function(p) {
    # paste0(
    # ifelse(ks.test(na.omit(p), y = "punif")$p.value < 0.05, "$\\star$ ", " "),

    formatC(100 * mean(p < 0.05, na.rm = TRUE), digits = 1, format = "f")
    # )
  }
)

mk <- simsalapar::mkAL(x = out_normal, vList = varList, repFirst = TRUE)
resArray <- getArray(mk, comp = "value")

# Compute percentage of rejection for a 5% test under H0
# Kolmogorov-Smirnov tests of uniformity
size_ar_normal <- apply(
  X = resArray[, 1, , , ],
  MARGIN = 1:3,
  FUN = function(p) {
    # paste0(
    # ifelse(ks.test(na.omit(p), y = "punif")$p.value < 0.05, "$\\star$ ", " "),

    formatC(100 * mean(p < 0.05, na.rm = TRUE), digits = 1, format = "f")
    # )
  }
)

size_ar <- array(
  dim = c(3, 3, 2, 3),
  dimnames = list(
    alt = c("$A_1$", "$A_2$", "$A_3$"),
    n = c("25", "50", "100"),
    distribution = c("Weibull", "normal"),
    m = c("2", "5", "10")
  )
)
size_ar[,, "Weibull", ] <- size_ar_weibull
size_ar[,, "normal", ] <- size_ar_normal
save(
  size_ar,
  file = "/home/lbelzile/Documents/Dropbox/Travail/Presentations/202604-Glasgow/results/sizeGEV.RData"
)

setwd(paste0(dirname(dirname(this.path::here())), "/results"))
write(
  print(
    toLatex(
      ftable(
        x = size_ar,
        col.vars = c("m", "n"),
        row.vars = c("distribution", "alt")
      ),
      ,
      booktabs = TRUE,
      center = TRUE,
      escapeFUN = identity,
      caption = "Size of tests (in percentage) at level 5\\% for the null hypothesis of max-stability when simulating from GEV penultimate approximations for different distributions and alternatives (alt). Stars indicate samples for which a Kolmogorov--Smirnov test of uniformity rejects the null hypothesis at level 5\\%."
    )
  ),
  file = "sizedistortion_GEV.tex",
  append = FALSE
)
