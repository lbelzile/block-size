library(mev)
library(ggplot2)
library(patchwork)

set.seed(202603)
theta <- 0.5
m <- seq(1, 20, length.out = 101)
# Simulate from Tavares first-order max-autoregressive process
# with unconditional standard Gumbel margins
y <- rmar1(n = 1000, theta = theta, shape = 0)
gg1_tav <- ggplot(
  data = data.frame(y = y, x = seq_along(y)),
  mapping = aes(y = y, x = x)
) +
  geom_point() +
  labs(x = "time", y = "") +
  theme_classic()

# Plot the extrapolation of the location
# based on max-stability, exact and accounting for extremal index
gg2_tav <- ggplot(
  data = data.frame(
    location = c(
      log1p((m - 1) * theta),
      log(m),
      log(m * 0.5)
    ),
    m = rep(m, length.out = 3 * length(m)),
    estimator = rep(
      x = c("true", "iid", "stationary"),
      each = length(m)
    )
  ),
  mapping = aes(
    x = m,
    y = location,
    col = estimator
  )
) +
  geom_line() +
  MetBrewer::scale_color_met_d(name = "Hiroshige") +
  theme_classic()
gg1_tav + gg2_tav
