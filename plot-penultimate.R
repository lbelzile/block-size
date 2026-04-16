library(mev)
library(ggplot2)
library(patchwork)

# Compute penultimate approximation for standard Gaussian
# over different block size m
m_seq <- seq(from = 30L, to = 30 * 15, by = 15L)
penult <- penultimate(
  family = 'norm',
  ddensF = function(x) {
    -x * dnorm(x)
  },
  m = m_seq,
  method = "bm",
  returnList = FALSE
)

# Compare with max-stability extrapolation from m=30
pars <- as.data.frame(
  t(sapply(seq_along(m_seq), function(m) {
    c(
      maxstable(
        pars = c(penult$loc[1], penult$scale[1], penult$shape[1]),
        m = m
      ),
      m = m_seq[m]
    )
  }))
)
# Build data frame for ggplot
penult_df <- data.frame(
  rbind(penult, pars)
)
penult_df$approximation <- factor(rep(
  c("penultimate", "max-stability"),
  each = length(m_seq)
))
g1a <- ggplot(
  data = penult_df,
  mapping = aes(x = m, y = loc, linetype = approximation)
) +
  geom_line() +
  labs(x = "block size", y = "location", subtitle = "normal")
g2a <- ggplot(
  data = penult_df,
  mapping = aes(x = m, y = scale, linetype = approximation)
) +
  geom_line() +
  labs(x = "block size", y = "scale")
g3a <- ggplot(
  data = penult_df,
  mapping = aes(x = m, y = shape, linetype = approximation)
) +
  geom_line() +
  labs(x = "block size", y = "shape")

# Same, this time with Weibull distribution
penult <- penultimate(
  family = 'weibull',
  shape = 0.8,
  m = m_seq,
  method = "bm",
  returnList = FALSE
)
pars <- as.data.frame(
  t(sapply(seq_along(m_seq), function(m) {
    c(
      maxstable(
        pars = c(
          penult$loc[1],
          penult$scale[1],
          penult$shape[1]
        ),
        m = m
      ),
      m = m_seq[m]
    )
  }))
)
penult_df <- data.frame(rbind(penult, pars))
penult_df$approximation <- factor(rep(
  c("penultimate", "max-stability"),
  each = length(m_seq)
))
g1b <- ggplot(
  data = penult_df,
  mapping = aes(x = m, y = loc, linetype = approximation)
) +
  geom_line() +
  labs(x = "block size", y = "location", subtitle = "Weibull")
g2b <- ggplot(
  data = penult_df,
  mapping = aes(x = m, y = scale, linetype = approximation)
) +
  geom_line() +
  labs(x = "block size", y = "scale")
g2c <- ggplot(
  data = penult_df,
  mapping = aes(x = m, y = shape, linetype = approximation)
) +
  geom_line() +
  labs(x = "block size", y = "shape")

# Pring the plots
(g1a + g2a + g3a) /
  (g1b + g2b + g2c) +
  plot_layout(guides = "collect") &
  theme_classic() + theme(legend.position = "bottom")
