# Simulation study for power
# of tests for block size

setwd(this.path::here())
source("helpers.R")
library(mev)
library(simsalapar)

date <- 20250505
seed_init <- floor(runif(1)[1] * date)

# Define some constants and settings
nobs_seq <- c(25, 50, 100)
B <- 2000L
m_seq <- c(2L, 5L, 10L)
id_seq <- c("weibull","normal")
icens_seq <- c(FALSE, TRUE)
lcens_seq <- c(0, 0.2)
# Cluster
ncores <- 50L
block <- 10L
delta_seq <- seq(1, 18, by = 0.25)

#Set list of variables for the simulation study
varList <- simsalapar::varlist(
  n.sim = list(type = "N", expr = quote(N), value = B),
  delta = list(type = "grid", value = delta_seq),
  nobs = list(type = "grid", value = nobs_seq),
  m = list(type = "grid", value = m_seq),
  id = list(type = "grid", value = id_seq),
  icens = list(type = "grid", value = icens_seq),
  lcens = list(type = "grid", value = lcens_seq),
  alt = list(type = "inner", value = c(1:3))
)

out <- doClusterApply(
  vList = varList,
  doAL = FALSE,
  sfile = paste0("power-study-blocksize-5.rds"),
  cluster = parallel::makeCluster(ncores, type = "PSOCK"),
  block.size = block,
  doOne = simu_fn_cens,
  check = FALSE,
  keepSeed = FALSE,
  seed = seed_init + 1:B,
  exports = ls()
)
mk <- simsalapar::mkAL(x = out, vList = varList, repFirst = TRUE)
power <- simsalapar::array2df(getArray(mk, comp = "value")) |>
  dplyr::group_by(alt, delta, nobs, m, id, icens, lcens) |>
  dplyr::summarize(
    power = mean(value < 0.05, na.rm = TRUE),
    .groups = "drop_last"
  ) |>
  dplyr::mutate(delta = delta_seq[as.integer(delta)])

save(power, file = paste0("power-study-5-power.RData"))
