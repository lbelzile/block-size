setwd(this.path::here())

library(mev)
library(lubridate)
library(xts)
library(ggplot2)
library(patchwork)

source("helpers.R")

data(abisko, package = "mev")

# Fill in series to replace missing days with zeros
# to do this, we create another time series of zeros
# and sum it with the original
# so as to compute 3 day sum of full time series
alldays <- seq.Date(abisko$date[1], to = tail(abisko$date, 1), by = 1)
zeros <- xts(x = rep(0, length(alldays)), order.by = alldays)
ts <- xts(x = abisko$precip, order.by = abisko$date)
cumul <- apply.daily(rbind(ts, zeros), sum)
names(cumul) <- "rain3day"
# Use rolling function to compute 3 day sum
abisko3 <- round(zoo::rollsum(cumul, 3), 3)

plot(
  abisko3,
  main = "Cumulative 3-day rainfall in Abisko (mm)",
  main.timespan = FALSE,
  format.labels = "%Y"
)

# Keep months June-August, and exactly 12 weeks a year
# The aggregate to four week period (m=4)
# And from there to yearly (m=3)
days_seq <- seq(
  from = yday(ymd("2026-06-15")),
  by = 1,
  length.out = 7 * 3 * 4
)
# Data series
xdat <- abisko3[yday(abisko3) %in% days_seq]
plot(
  x = (year(xdat) + yday(xdat) / diff(range(days_seq)))[xdat > 10],
  y = xdat[xdat > 10],
  ylab = "cumulative three-day rainfall",
  xlab = "day of year",
  bty = "l",
  ylim = c(10, max(xdat) + 1),
  yaxs = "i",
  pch = 20
)
# Make a data frame with consecutive time stamps
df_abisko3 <- data.frame(
  time = year(xdat) + yday(xdat) / diff(range(days_seq)),
  rain = as.numeric(xdat)
)

xdat_week <- build.blocks(
  df_abisko3$rain,
  block = 7,
  m = 4
)
xdat_month <- build.blocks(
  df_abisko3$rain,
  block = 7 * 4,
  m = 3
)
# Left-censoring bound
lb <- 10
mean(xdat_week < lb)
mean(xdat_month < lb)
mean(xdat_month[, 3] < lb)

# Estimate parameters with different block sizes
fit_week <- fit.gevblock(
  xdat = xdat_week,
  lb = lb,
  rounding = 0.1
)
fit_month <- fit.gevblock(
  xdat = xdat_month,
  lb = lb,
  rounding = 0.1
)
fit_year <- fit.gevblock(
  xdat = xdat_month[, 3, drop = FALSE],
  lb = lb,
  rounding = 0.1
)

# Likelihood ratio tests
test_week <- test.blocksize(
  xdat = xdat_week,
  lb = lb,
  rounding = 0.1
)
test_month <- test.blocksize(
  xdat = xdat_month,
  lb = lb,
  rounding = 0.1
)

# Quantile-quantile plots for the model

set.seed(202604)
qq_week <- qqplot.blocksize(
  xdat = xdat_week,
  lb = lb,
  type = c("all", "max"),
  B = 200,
  marginal = FALSE,
  rounding = 0.1,
  plot = FALSE,
  n = 250
)


qq_month <- qqplot.blocksize(
  xdat = xdat_month,
  lb = lb,
  type = c("all", "max"),
  B = 200,
  marginal = FALSE,
  rounding = 0.1,
  plot = FALSE,
  np = 250
)
autoplot(qq_week, type = "all") + autoplot(qq_month, type = "all")


# Check coefficients with max-stability
rbind(month2week = maxstable(fit_week, 4), month = fit_month)

# Profile likelihood

g1_abi <- ggplot(
  data = df_abisko3,
  mapping = aes(x = time, y = rain)
) +
  geom_point() +
  scale_y_continuous(limits = c(10, NA), expand = expansion()) +
  scale_x_continuous(
    breaks = c(1913, seq(1920, 2010, by = 10L), 2014),
    limits = c(1913, 2015),
    expand = expansion()
  ) +
  labs(x = "", y = "cumulative three-day rainfall (mm)") +
  theme_classic()

# Function to calculate the GEV log likelihood with censored data
gev.loglik.cens <- function(pars, xdat, lb, interval) {
  delta <- interval[1] / 2
  xdat <- as.numeric(xdat)
  if (pars[2] < 0 | pars[3] < -1) {
    return(-Inf)
  }
  obj <- sum(ifelse(
    xdat < lb,
    mev::pgev(lb, pars[1], pars[2], pars[3], log.p = TRUE),
    log(
      mev::pgev(xdat + delta, pars[1], pars[2], pars[3]) -
        mev::pgev(xdat - delta, pars[1], pars[2], pars[3])
    )
  ))
  obj
}

#' Profile log-likelihood for probability of exceedance
gev.prof.prob <- function(
  pars,
  psi,
  xquant,
  xdat,
  lb,
  interval,
  m = 1
) {
  sigma <- exp(pars[1])
  xi <- pars[2]
  mu <- xquant - sigma / xi * ((-log(1 - psi) / m)^(-xi) - 1)
  gev.loglik.cens(
    c(mu, sigma, xi),
    xdat = xdat,
    lb = lb,
    interval = interval
  )
}

# Sequence of values for log-probability at which to profile
psi_seq <- exp(seq(-8.5, -2, by = 0.05))
# Matrices for storing the other parameter vectors and the
# negative profile log likelihood values
prof_pars_mth <- matrix(ncol = 2, nrow = length(psi_seq))
prof_nll_mth <- numeric(length = length(psi_seq))

prof_pars_yr <- matrix(ncol = 2, nrow = length(psi_seq))
prof_nll_yr <- numeric(length = length(psi_seq))

xdat_year <- xdat_month[, 3]

# Obtain initial values for optimization from MLE
init_mth <- c(log(fit_month[2]), fit_month[3])
init_yr <- c(log(fit_year[2]), fit_year[3])
# Max-stability extrapolation of parameters from fit to monthly max to yearly
m2y <- maxstable(fit_month, m = 3)
# Maximum likelihood estimates from the models for annual probability of exceedance of 69.9 mm
pexc_mle_mth <- mev::pgev(
  q = 69.9,
  m2y[1],
  m2y[2],
  m2y[3],
  lower.tail = FALSE
)
pexc_mle_yr <- mev::pgev(
  q = 69.9,
  fit_year[1],
  fit_year[2],
  fit_year[3],
  lower.tail = FALSE
)

# Profile loop: start from MLE index, and go up/down
is <- which.min(abs(psi_seq - pexc_mle_mth))
for (i in is:length(psi_seq)) {
  opt_mth <- optim(
    par = init_mth,
    fn = gev.prof.prob,
    method = "N",
    control = list(fnscale = -1),
    psi = psi_seq[i],
    xquant = 69.9,
    xdat = xdat_month,
    lb = lb,
    m = 3,
    interval = 0.1
  )
  prof_nll_mth[i] <- opt_mth$value
  prof_pars_mth[i, ] <- opt_mth$par
  init_mth <- opt_mth$par
  opt_yr <- optim(
    par = init_yr,
    fn = gev.prof.prob,
    method = "N",
    control = list(fnscale = -1),
    psi = psi_seq[i],
    xquant = 69.9,
    xdat = xdat_month[, 3],
    lb = lb,
    m = 1,
    interval = 0.1
  )
  prof_nll_yr[i] <- opt_yr$value
  prof_pars_yr[i, ] <- opt_yr$par
  init_yr <- opt_yr$par
}
for (i in (is - 1):1) {
  opt_mth <- optim(
    par = init_mth,
    fn = gev.prof.prob,
    method = "N",
    control = list(fnscale = -1),
    psi = psi_seq[i],
    xquant = 69.9,
    xdat = xdat_month,
    lb = lb,
    m = 3,
    interval = 0.1
  )
  prof_nll_mth[i] <- opt_mth$value
  prof_pars_mth[i, ] <- opt_mth$par
  init_mth <- opt_mth$par
  opt_yr <- optim(
    par = init_yr,
    fn = gev.prof.prob,
    method = "N",
    control = list(fnscale = -1),
    psi = psi_seq[i],
    xquant = 69.9,
    xdat = xdat_month[, 3],
    lb = lb,
    m = 1,
    interval = 0.1
  )
  prof_nll_yr[i] <- opt_yr$value
  prof_pars_yr[i, ] <- opt_yr$par
  init_yr <- opt_yr$par
}

# Calculate the profile log likelihood at the MLE
max_pll_mth <- gev.loglik.cens(
  pars = fit_month,
  xdat = xdat_month,
  lb = lb,
  interval = 0.1
)
max_pll_yr <- gev.loglik.cens(
  pars = fit_year,
  xdat = xdat_year,
  lb = lb,
  interval = 0.1
)

# Create a list wrapper, to pass to 'mev' eprof
# routine to calculate confidence intervals
pll_mth <- list(
  pll = -prof_nll_mth,
  mle = pexc_mle_mth,
  maxpll = max_pll_mth,
  psi = psi_seq,
  psi.max = pexc_mle_mth,
  r = sign(pexc_mle_mth - psi_seq) *
    sqrt(2 * (max(prof_nll_mth) - prof_nll_mth))
)

pll_yr <- list(
  pll = -prof_nll_yr,
  mle = pexc_mle_yr,
  maxpll = max_pll_yr,
  psi = psi_seq,
  psi.max = pexc_mle_yr,
  r = sign(pexc_mle_yr - psi_seq) *
    sqrt(2 * (max(prof_nll_yr) - prof_nll_yr))
)
class(pll_mth) <- class(pll_yr) <- "eprof"
(confint(pll_mth))
(confint(pll_yr))
plot(
  psi_seq,
  prof_nll_mth - max_pll_mth,
  type = "l",
  xlab = "p",
  ylab = "profile log likelihood",
  ylim = c(-4, 0),
  xlim = c(0, 0.05),
  xaxs = "i",
  lwd = 1.5,
  panel.first = {
    abline(
      h = -qchisq(c(0.95, 0.99), df = 1) / 2,
      lty = 2,
      col = "grey"
    )
  },
  bty = "n"
)
lines(psi_seq, prof_nll_yr - max_pll_yr, col = 2, lwd = 1.5)
ns <- length(psi_seq)
conf_mth <- confint(pll_mth)
rug(conf_mth, lwd = 2)
conf_yr <- confint(pll_yr)
rug(conf_yr, col = 2, lwd = 2)

# Do the same plot, with ggplot2 instead
prof_abisko_df <-
  data.frame(
    p = rep(psi_seq, length.out = 2 * ns),
    profile = c(prof_nll_mth - max_pll_mth, prof_nll_yr - max_pll_yr),
    type = factor(rep(c("28-days", "yearly"), each = ns))
  )

g2_abi <- ggplot(
  data = prof_abisko_df,
  mapping = aes(x = p, col = type)
) +
  geom_hline(
    yintercept = -qchisq(c(0.95, 0.99), df = 1) / 2,
    alpha = 0.5,
    linetype = "dashed"
  ) +
  geom_line(mapping = aes(y = profile), linewidth = 1.25) +
  geom_rug(
    data = data.frame(
      p = c(conf_mth, conf_yr),
      type = factor(rep(c("28-days", "yearly"), each = 3))
    ),
    mapping = aes(x = p, col = type),
    linewidth = 1.25,
    sides = "b"
  ) +
  # coord_cartesian(clip = "off") +
  scale_x_continuous(
    limits = c(0, 0.05),
    expand = expansion(),
    oob = scales::oob_keep
  ) +
  scale_y_continuous(
    limits = c(-4, 0.1),
    expand = expansion(add = c(0, 0.1)),
    oob = scales::oob_keep
  ) +
  MetBrewer::scale_color_met_d("Hiroshige") +
  theme_classic()

g1_abi + g2_abi


# Bayesian analysis with default non informative priors
# not reported in the paper, yields similar results
gev.lpost <- function(pars, xdat, lb, interval) {
  revdbayes::gev_flat(pars, min_xi = -1) +
    gev.loglik.cens(
      pars = pars,
      xdat = xdat,
      lb = lb,
      interval = interval
    )
}

rou <- rust::ru(
  logf = gev.lpost,
  n = 10000,
  d = 3,
  init = fit_month,
  mode = fit_month,
  xdat = c(xdat_month),
  lb = lb,
  interval = 0.1,
  lower = c(-Inf, 1e-8, -1)
)$sim_vals

post <- apply(rou, 1, function(pars) {
  ypars <- maxstable(pars, 3)
  mev::pgev(q = 69.9, ypars[1], ypars[2], ypars[3], lower.tail = FALSE)
})
qupost <- quantile(post, probs = c(0.025, 0.5, 0.975))
qupost
