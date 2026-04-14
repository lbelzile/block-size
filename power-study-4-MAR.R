# Simulation study for power of tests for block size
# with maximum of blocks of size b from max-autoregressive processes

setwd(this.path::here())
source("helpers.R")
library(mev)
library(simsalapar)
date <- 20250414
seed_init <- floor(runif(4)[4] * date)
# Define some constants and settings
B <- 2000L
nobs_seq <- c(25L, 50L)
shape_seq <- c(0, 0.2, 0.4)
theta_seq <- c(seq(1, 0.25, by = -0.05))
block_seq <- c(2L, 4L, 8L)
# Cluster
ncores <- 50L
block <- 10L


# Run simulations
#Set list of variables for the simulation study
varList <- simsalapar::varlist(
  n.sim = list(type = "N", expr = quote(N), value = B),
  theta = list(type = "grid", value = theta_seq),
  nobs = list(type = "grid", value = nobs_seq),
  shape = list(type = "grid", value = shape_seq),
  block = list(type = "grid", value = block_seq),
  alt = list(type = "inner", value = 1:3),
  m = list(type = "frozen", value = 2L)
)

doClusterApply(
  vList = varList,
  doAL = FALSE,
  sfile = "power-study-blocksize-4.rds",
  cluster = parallel::makeCluster(ncores, type = "PSOCK"),
  block.size = block,
  doOne = simu_fn_st,
  keepSeed = FALSE,
  seed = seed_init + 1:B,
  exports = ls()
)
