setwd(this.path::here())
# Analysis of river flow of the Thames at Kingston
library(mev)
library(ggplot2)
library(exdex)

data(thames, package = "mev")
# Fit GEV model to m=1, 2 years max
flow <- thames$flow
flow2 <- build.blocks(
  thames$flow,
  block = 2,
  m = 1
)

# Parameter estimation via maximum likelihood
fit_m1 <- mev::fit.gev(flow)
fit_m1_gumbel <- mev::fit.gev(flow, fpar = list(shape = 0))
# Check for difference in fit
anova(fit_m1, fit_m1_gumbel)
# Check max-stability
(coef_ms2 <- mev::maxstable(coef(fit_m1), m = 2))
fit_m2 <- mev::fit.gev(flow2)
(coef_m2 <- coef(fit_m2))

# P-values for the test of max-stability with m=2
pvals_thames <- test.blocksize(
  xdat = build.blocks(
    thames$flow,
    block = 1,
    m = 2
  )
)
pvals_thames
# Check for extremal clustering
exdex::spm(data = flow, b = 3)

# Diagnostic plots
set.seed(2026)
qqplot.blocksize(
  xdat = build.blocks(
    thames$flow,
    block = 1,
    m = 2
  ),
  B = 200
)
