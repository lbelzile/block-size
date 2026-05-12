setwd(this.path::here())

library(mev)
library(lubridate)
library(xts)
library(ggplot2)
library(patchwork)

source("helpers.R")

# Profile likelihood

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

# Keep months June-September, and exactly 12 weeks a year
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

plot(
  x = yday(xdat),
  y = xdat,
  ylab = "cumulative three-day rainfall",
  xlab = "day of year",
  bty = "l",
  ylim = c(0, max(xdat) + 1),
  yaxs = "i",
  pch = 20,
  col = ifelse(xdat < 10, "grey", "black"),
  panel.last = {
    abline(h = 10, col = "grey")
  }
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

xdat_year <- xdat_month[, 3]

# Left-censoring bound
lb <- 5
mean(xdat_week < lb)
mean(xdat_month < lb)
mean(xdat_month[, 3] < lb)

# Preliminary optimization - without rounding
fit_week_prelim <- fit.gevblock(xdat = xdat_week, lb = lb)
fit_month_prelim <- fit.gevblock(xdat = xdat_month, lb = lb)
fit_year_prelim <- fit.gevblock(
  xdat = xdat_month[, 3, drop = FALSE],
  lb = lb
)

# Estimate parameters with different block sizes
fit_week <- optim(
  fn = gev.loglik.cens,
  par = fit_week_prelim,
  xdat = xdat_week,
  lb = lb,
  rounding = 0.1,
  control = list(fnscale = -1),
  method = "Nelder")$par
fit_month <- optim(
  fn = gev.loglik.cens,
  par = fit_month_prelim,
  xdat = xdat_month,
  lb = lb,
  rounding = 0.1,
  control = list(fnscale = -1),
  method = "Nelder")$par
fit_year <- optim(
  fn = gev.loglik.cens,
  par = fit_year_prelim,
  xdat = xdat_year,
  lb = lb,
  rounding = 0.1,
  control = list(fnscale = -1),
  method = "Nelder")$par

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

# Sensitivity analysis

test.blocksize(xdat = xdat_week, lb = 5, rounding = 0.1)
test.blocksize(xdat = xdat_week, lb = 15, rounding = 0.1)

test.blocksize(xdat = xdat_month, lb = 5, rounding = 0.1)
test.blocksize(xdat = xdat_month, lb = 15, rounding = 0.1)

# Compare parameter estimates when modifying the threshold
optim(
  fn = gev.loglik.cens,
  par = fit_week_prelim,
  xdat = xdat_week,
  lb = 5,
  rounding = 0.1,
  control = list(fnscale = -1),
  method = "Nelder")$par

optim(
  fn = gev.loglik.cens,
  par = fit_week_prelim,
  xdat = xdat_week,
  lb = 15,
  rounding = 0.1,
  control = list(fnscale = -1),
  method = "Nelder")$par

# Quantile-quantile plots for the model

set.seed(202604)
qq_week <- qqplot.blocksize(
  xdat = xdat_week,
  lb = lb,
  type = c("all", "max"),
  B = 2000,
  marginal = FALSE,
  rounding = 0.1,
  plot = FALSE,
  n = 250
)

set.seed(202604)
qq_month <- qqplot.blocksize(
  xdat = xdat_month,
  lb = lb,
  type = c("all", "max"),
  B = 2000,
  marginal = FALSE,
  rounding = 0.1,
  plot = FALSE,
  np = 250
)
save(qq_week, qq_month, file = "results/qqplots-Abisko.RData")
autoplot(qq_week, type = "all") +
  labs(caption = "7-days (weekly) maximum") +
  autoplot(qq_month, type = "all") +
  labs(caption = "28-days (monthly) maximum")
ggsave(filename = "fig/Abisko-pp-plots.pdf", width = 10, height = 6, units = "in")

# Check coefficients with max-stability
rbind(month2week = maxstable(fit_week, 4), month = fit_month)

# Sequence of values for log-probability at which to profile
psi_seq <- exp(seq(-8.5, -2, by = 0.05))
# Matrices for storing the other parameter vectors and the
# negative profile log likelihood values
prof_pars_wk <- prof_pars_mth <- prof_pars_yr <- matrix(
  ncol = 2,
  nrow = length(psi_seq)
)
prof_nll_wk <- prof_nll_mth <- prof_nll_yr <- numeric(length = length(psi_seq))


# Obtain initial values for optimization from MLE
init_wk <- c(log(fit_week[2]), fit_week[3])
init_mth <- c(log(fit_month[2]), fit_month[3])
init_yr <- c(log(fit_year[2]), fit_year[3])
# Max-stability extrapolation of parameters from fit to monthly max to yearly
w2y <- maxstable(fit_week, m = 12)
m2y <- maxstable(fit_month, m = 3)
# Maximum likelihood estimates from the models for annual probability of exceedance of 69.9 mm

pexc_mle_wk <- mev::pgev(
  q = 69.9,
  w2y[1],
  w2y[2],
  w2y[3],
  lower.tail = FALSE
)
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
n_psi <- length(psi_seq)
for (i in c(is:n_psi, (is - 1):1)) {
  opt_wk <- optim(
    par = init_wk,
    fn = gev.prof.prob,
    method = "N",
    control = list(fnscale = -1),
    psi = psi_seq[i],
    xquant = 69.9,
    xdat = xdat_week,
    lb = lb,
    m = 12,
    rounding = 0.1
  )
  prof_nll_wk[i] <- opt_wk$value
  prof_pars_wk[i, ] <- opt_wk$par
  if (i != n_psi) {
    init_wk <- opt_wk$par
  } else {
    init_wk <- prof_pars_wk[is, ]
  }
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
    rounding = 0.1
  )
  prof_nll_mth[i] <- opt_mth$value
  prof_pars_mth[i, ] <- opt_mth$par
  if (i != n_psi) {
    init_mth <- opt_mth$par
  } else {
    init_mth <- prof_pars_mth[is, ]
  }
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
    rounding = 0.1
  )
  prof_nll_yr[i] <- opt_yr$value
  prof_pars_yr[i, ] <- opt_yr$par
  if (i != n_psi) {
    init_yr <- opt_yr$par
  } else {
    init_yr <- prof_pars_yr[is, ]
  }
}
# Calculate the profile log likelihood at the MLE
max_pll_wk <- gev.loglik.cens(
  pars = fit_week,
  xdat = xdat_week,
  lb = lb,
  rounding = 0.1
)
max_pll_mth <- gev.loglik.cens(
  pars = fit_month,
  xdat = xdat_month,
  lb = lb,
  rounding = 0.1
)
max_pll_yr <- gev.loglik.cens(
  pars = fit_year,
  xdat = xdat_year,
  lb = lb,
  rounding = 0.1
)

# Create a list wrapper, to pass to 'mev' eprof
# routine to calculate confidence intervals
pll_wk <- list(
  pll = -prof_nll_wk,
  mle = pexc_mle_wk,
  maxpll = max_pll_wk,
  psi = psi_seq,
  psi.max = pexc_mle_wk,
  r = sign(pexc_mle_wk - psi_seq) *
    sqrt(2 * (max(prof_nll_wk) - prof_nll_wk))
)

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
class(pll_wk) <- class(pll_mth) <- class(pll_yr) <- "eprof"
(confint(pll_wk))
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
lines(psi_seq, prof_nll_wk - max_pll_wk, col = 4, lwd = 1.5)

conf_wk <- confint(pll_wk)
rug(conf_wk, lwd = 2, col = 4)
conf_mth <- confint(pll_mth)
rug(conf_mth, lwd = 2)
conf_yr <- confint(pll_yr)
rug(conf_yr, col = 2, lwd = 2)

# Do the same plot, with ggplot2 instead
prof_abisko_df <-
  data.frame(
    p = rep(psi_seq, length.out = 3 * n_psi),
    profile = c(
      prof_nll_wk - max_pll_wk,
      prof_nll_mth - max_pll_mth,
      prof_nll_yr - max_pll_yr
    ),
    type =
      factor(
        rep(c("weekly", "monthly", "seasonal"), each = n_psi),
      levels = c("weekly", "monthly", "seasonal")
    )
  )


g1_abi <- ggplot(
  data = df_abisko3 |> dplyr::filter(rain > 0),
  mapping = aes(
    x = time,
    y = rain,
    col = factor(rain <= 10)
  )
) +
  geom_hline(yintercept = 10, col = "grey") +
  geom_point() +
  scale_color_grey() +
  scale_y_continuous(limits = c(0, NA), expand = expansion()) +
  scale_x_continuous(
    breaks = seq(1920, 2010, by = 10L),
    limits = c(1913, 2016),
    expand = expansion()
  ) +
  labs(
    x = "",
    y = "",
    subtitle = "cumulative three-day rainfall (mm)"
  ) +
  theme_classic() +
  theme(legend.position = "grey")


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
      p = c(conf_wk, conf_mth, conf_yr),
      type = factor(rep(c("weekly", "monthly", "seasonal"), each = 3),
                    levels = c("weekly", "monthly", "seasonal"))
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
  theme_classic() +
  theme(legend.position = "bottom")

g1_abi + g2_abi
ggsave(filename = "fig/Abisko-data-profile.pdf", width = 11, height = 6, units = "in")

# Bayesian analysis with default non informative priors
# not reported in the paper, yields similar results
gev.lpost <- function(pars, xdat, lb, rounding) {
  revdbayes::gev_flat(pars, min_xi = -1) +
    gev.loglik.cens(
      pars = pars,
      xdat = xdat,
      lb = lb,
      rounding = rounding
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
  rounding = 0.1,
  lower = c(-Inf, 1e-8, -1)
)$sim_vals

post <- apply(rou, 1, function(pars) {
  ypars <- maxstable(pars, 3)
  mev::pgev(q = 69.9, ypars[1], ypars[2], ypars[3], lower.tail = FALSE)
})
qupost <- quantile(post, probs = c(0.025, 0.5, 0.975))
qupost
