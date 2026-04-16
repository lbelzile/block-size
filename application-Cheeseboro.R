setwd(this.path::here())

library(dplyr)
library(mev)
library(ggplot2)
library(exdex)
library(lubridate)
library(patchwork)
library(xts)
library(tinytable)

source("helpers.R")

## Application 2
data(cheeseborowind, package = "mev")
# Change units from m/s to miles per hour (mph)
# since instrumental precision is +/- 1 mph
gust <- cheeseborowind$gust <- round(2.236936 * cheeseborowind$gust, 1)
gust2 <- build.blocks(gust, block = 2, m = 4)
gust4 <- build.blocks(gust, block = 4, m = 4)

# Calculate yearly maximum for plot
yearly_max <- xts::apply.yearly(
  xts::xts(gust, order.by = cheeseborowind$date),
  max
)

# Time series plot of wind gust speed as a function of day of the year
# plus rugs for the yearly maxima
gg_ts_ymax <- ggplot(
  data = cheeseborowind,
  mapping = aes(x = lubridate::yday(date), y = gust)
) +
  geom_point(position = "jitter") +
  geom_rug(
    data = as.data.frame(yearly_max),
    mapping = aes(x = NA, y = V1),
    sides = "r"
  ) +
  labs(
    x = "Day of year",
    y = "",
    subtitle = "Maximum daily wind gust speed (miles per hour)"
  ) +
  scale_x_continuous(limits = c(1, 31), expand = expansion(add = c(0.1, 1))) +
  theme_classic()

gg_ts_ymax

# P-values for the test of max-stability
# with different alternatives, block sizes and m
pvals_chd <- array(
  dim = c(4, 3, 2),
  dimnames = list(
    block = c("1", "2", "4", "8"),
    alternative = c("A1", "A2", "A3"),
    m = c("2", "4")
  )
)
pvals_chd[,, 1] <-
  rbind(
    test.blocksize(
      rounding = 1,
      xdat = build.blocks(gust, block = 1, m = 2)
    )$pval,
    test.blocksize(
      rounding = 1,
      xdat = build.blocks(gust, block = 2, m = 2)
    )$pval,
    test.blocksize(
      rounding = 1,
      xdat = build.blocks(gust, block = 4, m = 2)
    )$pval,
    test.blocksize(
      rounding = 1,
      xdat = build.blocks(gust, block = 8, m = 2)
    )$pval
  )
pvals_chd[,, 2] <- rbind(
  test.blocksize(
    rounding = 1,
    xdat = build.blocks(gust, block = 1, m = 4)
  )$pval,
  test.blocksize(
    rounding = 1,
    xdat = build.blocks(gust, block = 2, m = 4)
  )$pval,
  test.blocksize(
    rounding = 1,
    xdat = build.blocks(gust, block = 4, m = 4)
  )$pval,
  test.blocksize(
    rounding = 1,
    xdat = build.blocks(gust, block = 8, m = 4)
  )$pval
)

df_pvals_dch <- array2DF(pvals_chd, responseName = "pvalue") |>
  dplyr::mutate(
    block = factor(block),
    alternative = factor(alternative),
    m = factor(m)
  )


# Create LaTeX table for the manuscript
dfc_tab <- data.frame(
  block = c("\\(b=1\\)", "\\(b=2\\)", "\\(b=4\\)", "\\(b=8\\)"),
  round(cbind(pvals_chd[,, 1], pvals_chd[,, 2]), 3)
)
# dfc_tab
tt_dfc <- tt(dfc_tab)
names(tt_dfc) <- c(
  "block",
  "\\(A_1\\)",
  "\\(A_2\\)",
  "\\(A_3\\)",
  "\\(A_1\\)",
  "\\(A_2\\)",
  "\\(A_3\\)"
)
tt_dfc |>
  format_tt(j = 2:7, digits = 2) |>
  format_tt(escape = FALSE) |>
  group_tt(
    j = list(
      "\\(m=2\\)" = 2:4,
      "\\(m=4\\)" = 5:6
    )
  )


## Probability-probability plots

set.seed(202504)
qq2 <- qqplot.blocksize(
  xdat = gust2,
  type = c("max", "all"),
  marginal = TRUE,
  rounding = 1,
  B = 200,
  plot = FALSE,
  np = 125
)


set.seed(202604)
qq4 <- qqplot.blocksize(
  xdat = gust4,
  type = c("max", "all"),
  marginal = TRUE,
  rounding = 1,
  B = 200,
  plot = FALSE,
  np = 125
)

autoplot(qq2, type = "max") + autoplot(qq4, type = "max")

## Investigation of extremal dependence

xacd_b4 <- mev::xacf(
  gust,
  qlev = 0.95,
  B = 1000,
  confint = TRUE
)

# Selection of block size for Northrop (2015) estimator
b_select <- choose_b(
  data = gust,
  b = 2:16,
  bias_adjust = "N",
  interval_type = "lik",
  conf_scale = "theta"
)
plot(b_select)
# Extremal index
theta_b4_spm <- b_select$theta_sl[4, "N2015"]

theta_mle <- mev::xdep.xindex(
  xdat = gust,
  qlev = 0.95,
  estimator = "mle",
  confint = "lrt"
)

mle_b4 <- mev::fit.gevblock(xdat = gust4, rounding = 1)
med_b4s <- gev.Nyr(
  par = mle_b4,
  nobs = length(gust4),
  N = 50 * 32 / 4 * theta_b4_spm,
  type = "median"
)
med_b4 <- gev.Nyr(
  par = mle_b4,
  nobs = length(gust4),
  N = 50 * 32 / 4,
  type = "median"
)

# Median and mean of 50 year maximum
mle_med50_b4 <- gev.mle(
  xdat = c(gust4),
  args = c("Nmean", "Nquant"),
  N = 50 * 32 / 4,
  q = 0.5
)
# Profile likelihood estimation for the latter
prof <- gev.pll(
  psi = seq(80, 250, length.out = 101),
  param = "Nquant",
  dat = c(gust4),
  q = 0.5,
  N = 50 * 32 / 4,
  plot = FALSE
)
# 95% profile-based confidence intervals
conf_prof <- confint(prof)

# Gumbel plot of maxima of 1 year, 2 years, etc.
# Mimicking Cox, Isham and Northrop (2002)
n_ch <- length(gust)
df <- rbind(
  data.frame(x = -log(-log(ppoints(n_ch))), y = sort(gust), block = 1),
  data.frame(
    x = -log(-log(ppoints(n_ch / 2))),
    y = sort(build.blocks(gust, block = 2, m = 1)),
    block = 2
  ),
  data.frame(
    x = -log(-log(ppoints(n_ch / 4))),
    y = sort(build.blocks(gust, block = 4, m = 1)),
    block = 4
  ),
  data.frame(
    x = -log(-log(ppoints(n_ch / 8))),
    y = sort(build.blocks(gust, block = 8, m = 1)),
    block = 8
  ),
  data.frame(
    x = -log(-log(ppoints(n_ch / 16))),
    y = sort(build.blocks(gust, block = 16, m = 1)),
    block = 16
  ),
  data.frame(
    x = -log(-log(ppoints(n_ch / 32))),
    y = sort(build.blocks(gust, block = 32, m = 1)),
    block = 32
  )
)
df$block <- as.factor(df$block)

gg_gumbelplot <- ggplot() +
  geom_point(
    data = df,
    mapping = aes(x = x, y = y, col = block, shape = block)
  ) +
  MetBrewer::scale_color_met_d("Hiroshige") +
  labs(
    x = "theoretical Gumbel positions",
    subtitle = "daily maximum gust speed (mph)",
    y = "",
    colour = "block size",
    shape = "block size"
  ) +
  theme_classic() +
  theme(legend.position = "inside", legend.position.inside = c(0.9, 0.3))

# Comparing extrapolations of different
# fit to GEV to 50-year maxima
fit1 <- maxstable(
  fit.gevblock(build.blocks(gust), rounding = 1),
  m = 32 * 50
)
fit2 <- maxstable(
  fit.gevblock(build.blocks(gust, block = 2), rounding = 1),
  m = 16 * 50
)
fit4 <- maxstable(
  fit.gevblock(build.blocks(gust, block = 4), rounding = 1),
  m = 8 * 50
)
fit8 <- maxstable(
  fit.gevblock(build.blocks(gust, block = 8), rounding = 1),
  m = 4 * 50
)
fit16 <- maxstable(
  fit.gevblock(build.blocks(gust, block = 16), rounding = 1),
  m = 2 * 50
)


gg_50yrdens <- ggplot() +
  stat_function(
    fun = function(x) {
      mev::dgev(x, loc = fit1[1], scale = fit1[2], shape = fit1[3])
    },
    aes(col = "1"),
    linewidth = 1.5,
    n = 1001
  ) +
  stat_function(
    fun = function(x) {
      mev::dgev(x, loc = fit2[1], scale = fit2[2], shape = fit2[3])
    },
    aes(col = "2"),
    linewidth = 1.5,
    n = 1001
  ) +
  stat_function(
    fun = function(x) {
      mev::dgev(x, loc = fit4[1], scale = fit4[2], shape = fit4[3])
    },
    aes(col = "4"),
    linewidth = 1.5,
    n = 1001
  ) +
  stat_function(
    fun = function(x) {
      mev::dgev(x, loc = fit8[1], scale = fit8[2], shape = fit8[3])
    },
    aes(col = "8"),
    linewidth = 1.5,
    n = 1001
  ) +
  stat_function(
    fun = function(x) {
      mev::dgev(x, loc = fit16[1], scale = fit16[2], shape = fit16[3])
    },
    aes(col = "16"),
    linewidth = 1.5,
    n = 1001
  ) +
  scale_y_continuous(limits = c(0, 0.2), expand = expansion()) +
  scale_x_continuous(limits = c(100, 370) * 0.6, expand = expansion()) +
  scale_color_manual(
    name = 'block size (days)',
    breaks = paste0(c(1L, 2L, 4L, 8L, 16L)),
    labels = paste0(c(1L, 2L, 4L, 8L, 16L)),
    values = MetBrewer::met.brewer(name = "Hiroshige", n = 5)
  ) +
  labs(
    y = "",
    subtitle = "density of 50 year maximum",
    x = "wind gust speed (mph)"
  ) +
  theme_classic() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.8)
  )


gg_gumbelplot + gg_50yrdens
