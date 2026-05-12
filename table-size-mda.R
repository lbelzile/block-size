library(simsalapar)

B <- 2000L
delta_seq <- c(seq(1, 4, by = 0.5), seq(5, 16, by = 1))
nobs_seq <- c(25L, 50L, 100L)
id_seq <- 2:3
m_seq <- c(2L, 5L, 10L)

# Run simulations
#Set list of variables for the simulation study
varList <- simsalapar::varlist(
  n.sim = list(type = "N", expr = quote(N), value = B),
  delta = list(type = "grid", value = delta_seq),
  nobs = list(type = "grid", value = nobs_seq),
  id = list(type = "grid", value = id_seq),
  m = list(type = "grid", value = m_seq),
  m0 = list(type = "frozen", value = 30),
  alt = list(type = "inner", value = 1:3)
)
out <- readRDS("power-study-blocksize-2.rds")
mk <- simsalapar::mkAL(x = out, vList = varList, repFirst = TRUE)
resArray <- getArray(mk, comp = "value")

# Compute percentage of rejection for a 5% test under H0
# Kolmogorov-Smirnov tests of uniformity
size_ar <- apply(
  X = resArray[, "1.0", , , , ],
  MARGIN = 1:4,
  FUN = function(p) {
    paste0(
      ifelse(ks.test(na.omit(p), y = "punif")$p.value < 0.05, "$\\star$ ", " "),

      formatC(100 * mean(p < 0.05, na.rm = TRUE), digits = 1, format = "f")
    )
  }
)

dimnames(size_ar) <- list(
  alt = c("$A_1$", "$A_2$", "$A_3$"),
  n = c("25", "50", "100"),
  distribution = c("Weibull", "normal"),
  m = c("2", "5", "10")
)


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
      caption = "Size of tests (in percentage) at level 5\\% for the null hypothesis of max-stability when simulating from the max-domain of attraction for different distributions and alternatives (alt). Stars indicate samples for which a Kolmogorov--Smirnov test of uniformity rejects the null hypothesis at level 5\\%."
    )
  ),
  file = "sizedistortion_MDA.tex",
  append = FALSE
)
