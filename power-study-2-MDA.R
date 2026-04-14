# Simulation study for power of tests for block size
# with data drawn from max domain of attraction

setwd(this.path::here())
source("helpers.R")
library(mev)
library(simsalapar)
date <- 20250414
seed_init <- floor(runif(2)[2] * seed_init)

# Define some constants and settings
nobs_seq <- c(25L, 50L, 100L)
B <- 2000L
m_seq <- c(2L, 5L, 10L)
id_seq <- 2:3
delta_seq <- c(seq(1, 4, by = 0.5), seq(5, 16, by = 1))

# Cluster
ncores <- 50L
block <- 10L


# Run simulations
#Set list of variables for the simulation study
varList <- simsalapar::varlist(
  n.sim = list(type = "N", expr = quote(N), value = B),
  delta = list(type = "grid", value = delta_seq),
  nobs = list(type = "grid", value = nobs_seq),
  id = list(type = "grid", value = id_seq),
  m = list(type = "grid", value = m_seq),
  alt = list(type = "inner", value = 1:3),
  m0 = list(type = "frozen", value = 30)
)

doClusterApply(
  vList = varList,
  doAL = FALSE,
  sfile = "power-study-blocksize-2.rds",
  cluster = parallel::makeCluster(ncores, type = "PSOCK"),
  block.size = block,
  doOne = simu_fn_mda,
  keepSeed = FALSE,
  seed = seed_init + 1:B,
  exports = ls()
)
