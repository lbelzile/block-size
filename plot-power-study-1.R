# Simulation study for power
# of tests for block size

g1 <- ggplot(
  data = res1_df,
  mapping = aes(
    x = delta,
    y = 100 * power,
    col = alternative,
    linetype = nobs
  )
) +
  geom_hline(linetype = "dashed", alpha = 0.5, yintercept = 5) +
  geom_line() +
  labs(
    x = expression(delta),
    linetype = "n",
    y = "power (%)",
    subtitle = "Generalized extreme value distribution"
  ) +
  facet_wrap(facets = vars(m), labeller = label_both) +
  MetBrewer::scale_color_met_d(name = "Hiroshige") +
  theme_minimal() +
  theme(legend.position = 'bottom')

g2 <- ggplot(
  data = res2_df,
  mapping = aes(
    x = delta,
    y = 100 * power,
    col = alternative,
    linetype = nobs
  )
) +
  geom_hline(linetype = "dashed", alpha = 0.5, yintercept = 5) +
  geom_line() +
  facet_wrap(facets = vars(m), labeller = label_both) +
  labs(
    x = expression(delta),
    linetype = "n",
    y = "power (%)",
    subtitle = "GEV approximation to Weibull distribution"
  ) +
  MetBrewer::scale_color_met_d(name = "Hiroshige") +
  theme_minimal() +
  theme(legend.position = 'bottom')


g3 <- ggplot(
  data = res3_df,
  mapping = aes(
    x = delta,
    y = 100 * power,
    col = alternative,
    linetype = nobs
  )
) +
  geom_hline(linetype = "dashed", alpha = 0.5, yintercept = 5) +
  geom_line() +
  facet_wrap(facets = vars(m), labeller = label_both) +
  labs(
    x = expression(delta),
    linetype = "n",
    y = "power (%)",
    subtitle = "GEV approximation to normal distribution"
  ) +
  MetBrewer::scale_color_met_d(name = "Hiroshige") +
  theme_minimal() +
  theme(legend.position = 'bottom')

(g1 / g2 / g3) +
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')
ggsave(
  filename = "fig/power-curve-GEV-1.pdf",
  width = 10,
  height = 8,
  units = "in"
)
