library(mev)
library(ggplot2)
library(patchwork)
setwd(this.path::here())
# Comparison of the density of G^m (full black) versus the penultimate
# approximation for m=30 (dashed blue), The bottom panel shows the
# max-stable extrapolation of the latter (dashed blue) against the
# penultimate for m=300 (dotted orange).

dmax <- function(x, m, family = "norm", args = NULL) {
  exp(
    log(m) +
      (m - 1) *
      do.call(
        paste0("p", family),
        args = c(list(q = x, log.p = TRUE), args)
      ) +
      do.call(paste0("d", family), args = c(list(x = x, log = TRUE), args))
  )
}
m <- 30
xseq <- seq(0, 6, length.out = 1001)
cols <- c("grey10","grey50","grey80")
cols <- MetBrewer::met.brewer(name = "Hiroshige", n = 3)
dnorm_max <- dmax(
  x = xseq,
  m = m,
  family = "norm"
)

# Penultimate approximation to m=30 max
penult_norm <- penultimate(
  family = "norm",
  method = "bm",
  m = m,
  ddensF = function(x) {
    -x * dnorm(x)
  }
)
penult_weibull <- penultimate(
  family = "weibull",
  method = "bm",
  m = m,
  shape = 0.8
)

# Max-stable extrapolation and penultimate for m=300
e_norm300 <- maxstable(penult_norm[1:3], m = 10)
p_norm300 <- penultimate(
  family = "norm",
  method = "bm",
  m = 10*m,
  ddensF = function(x) {
    -x * dnorm(x)
  }
)

e_weib300 <- maxstable(penult_weibull[1:3], m = 10)
p_weib300 <- penultimate(
  family = "weibull",
  method = "bm",
  m = 10*m,
  shape = 0.8
)

# Create plots
xlim <- c(0.5, 5)
gnorm <- ggplot() +
  stat_function(
    fun = dmax,
    args = list(family = "norm", m = m),
    n = 1001,
    linewidth = 1.01,
    # alpha = 0.7,
    xlim = xlim
  ) +
  stat_function(
    fun = dgev,
    args = penult_norm[1:3],
    color = cols[3],
    linewidth = 1.01,
    n = 1001,
    # alpha = 0.7,
  ) +
  stat_function(
    fun = dmax,
    args = list(family = "norm", m = m*10),
    n = 1001,
    xlim = xlim,
    linetype = "longdash"
  ) +
  stat_function(
    fun = dgev,
    args = e_norm300[1:3],
    color = cols[3],
    n = 1001,
    linetype = "longdash"
  ) +
  # stat_function(
  #   fun = dgev,
  #   args = p_norm300[1:3],
  #   col = cols[2],
  #   n = 1001,
  #   linewidth = 1.1,
  #   linetype = "longdash"
  # ) +
  scale_x_continuous(
    limits = xlim,
    expand = expansion()
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(add = c(0, 0.01))
  ) +
  labs(y = "density", caption = "normal") +
  theme_classic()
xlim <- c(0, 20)
gweibull <- ggplot() +
  stat_function(
    fun = dmax,
    args = list(
      family = "weibull",
      m = m,
      args = list(shape = 0.8)
    ),
    linewidth = 1.01,
    n = 1001,
    xlim = xlim
  ) +
  stat_function(
    fun = dgev,
    args = penult_weibull[1:3],
    col = cols[3],
    n = 1001,
    # linetype = "longdash",
    linewidth = 1.01
  ) +
  stat_function(
    fun = dmax,
    args = list(family = "weibull", m = m*10, args = list(shape = 0.8)),
    n = 1001,
    linetype = "longdash",
    # linewidth = 1.01,
    xlim = xlim
  ) +
  stat_function(
    fun = dgev,
    args = e_weib300[1:3],
    col = cols[3],
    n = 1001,
    # linewidth = 1.01,
    linetype = "longdash"
  ) +
  # stat_function(
  #   fun = dgev,
  #   args = p_weib300[1:3],
  #   col = cols[2],
  #   n = 1001,
  #   linewidth = 1.1,
  #   linetype = "longdash"
  # ) +
  scale_x_continuous(
    limits = xlim,
    expand = expansion()
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(add = c(0, 0.01))
  ) +
  labs(y = "density", caption = "Weibull(0.8)") +
  theme_classic()

gweibull + gnorm

ggsave("fig/penultimate-density.pdf",
       width = 10, height = 4, units = "in")

