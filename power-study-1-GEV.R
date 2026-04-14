# Simulation study for power
# of tests for block size

setwd(this.path::here())
source("helpers.R")
library(mev)
library(simsalapar)

date <- 20250414
seed_init <- floor(runif(1)[1] * date)

# Define some constants and settings
nobs_seq <- c(25, 50, 100)
B <- 2000L
m_seq <- c(2L, 5L, 10L)
modelid <- 1:3
# Cluster
ncores <- 50L
block <- 10L

# Run simulations
for (id in modelid) {
  if (id == 1L) {
    delta_seq <- seq(0, 4, by = 0.25)
  } else {
    delta_seq <- seq(1, 18, by = 0.25)
  }
  #Set list of variables for the simulation study
  varList <- simsalapar::varlist(
    n.sim = list(type = "N", expr = quote(N), value = B),
    delta = list(type = "grid", value = delta_seq),
    nobs = list(type = "grid", value = nobs_seq),
    m = list(type = "grid", value = m_seq),
    id = list(type = "frozen", value = id),
    alt = list(type = "inner", value = c(1:3))
  )

  doClusterApply(
    vList = varList,
    doAL = FALSE,
    sfile = paste0("power-study-blocksize-1_", id, ".rds"),
    cluster = parallel::makeCluster(ncores, type = "PSOCK"),
    block.size = block,
    doOne = simu_fn,
    keepSeed = FALSE,
    seed = seed_init + 1:B,
    exports = ls()
  )
}
