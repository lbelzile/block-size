# Simulation study for power of tests for block size
# from first-order max-autoregressive processes

setwd(this.path::here())
source("helpers.R")
library(mev)
library(simsalapar)
date <- 20250414
seed_init <- floor(runif(3)[3] * date)
# Define some constants and settings
nobs_seq <- c(25L, 50L, 100L)
B <- 2000L
m_seq <- c(2L, 5L, 10L)
shape_seq <- c(0, 0.2, 0.4)
theta_seq <- c(seq(1, 0.2, by = -0.05))

# Cluster
ncores <- 40L
block <- 10L


# Run simulations
#Set list of variables for the simulation study
varList <- simsalapar::varlist(
  n.sim = list(type = "N", expr = quote(N), value = B),
  theta = list(type = "grid", value = theta_seq),
  nobs = list(type = "grid", value = nobs_seq),
  shape = list(type = "grid", value = shape_seq),
  m = list(type = "grid", value = m_seq),
  alt = list(type = "inner", value = 1:3),
  block = list(type = "frozen", value = 1L)
)

out <- doClusterApply(
  vList = varList,
  doAL = FALSE,
  sfile = "power-study-blocksize-3.rds",
  cluster = parallel::makeCluster(ncores, type = "PSOCK"),
  block.size = block,
  doOne = simu_fn_st,
  check = FALSE,
  keepSeed = FALSE,
  seed = seed_init + 1:B,
  exports = ls()
)

mk <- simsalapar::mkAL(x = out, vList = varList, repFirst = TRUE)
power3 <- simsalapar::array2df(getArray(mk, comp = "value")) |>
  dplyr::group_by(alt, theta, shape, nobs, m) |>
  dplyr::summarize(power = mean(value < 0.05, na.rm = TRUE)) |>
  dplyr::mutate(theta = theta_seq[as.integer(theta)])
save(power3, file = paste0("power-study-3.RData"))
