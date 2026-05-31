# Simulation study for power
# of tests for block size
setwd(this.path::here())
load("results/power-study-5-power.RData")

library(dplyr)
library(ggplot2)
library(patchwork)



g1 <- ggplot(
  data = power |>
    dplyr::filter(id == "weibull",  m == 2, delta < 5
                  # (delta < 5 & m == 2) |
                  # (delta < 13 & m == 5) |
                  #   m == 10
                  ) |>
    dplyr::mutate(alternative = alt,
                  cens = factor(
      dplyr::case_when((icens == "FALSE" & lcens == "0.0") ~ "none",
                (icens == "FALSE" & lcens == "0.2") ~ "left-censored",
                (icens == "TRUE" & lcens == "0.0") ~ "rounded",
                (icens == "TRUE" & lcens == "0.2") ~ "censored"))
      ),
  mapping = aes(
    x = delta,
    y = 100 * power,
    linetype = nobs,
    col = cens
  )
) +
  geom_hline(linetype = "dashed", alpha = 0.5, yintercept = 5) +
  geom_line() +
  labs(
    x = expression(delta),
    linetype = "sample size",
    color = "censoring",
    y = "power (%)",
    subtitle = "Generalized extreme value approximation to Weibull(0.8)"
  ) +
  facet_wrap(
    facets = vars(alternative),
    labeller = label_both,
    scales = "free_x") +
  MetBrewer::scale_color_met_d(name = "Hiroshige") +
  theme_minimal() +
  theme(legend.position = 'bottom')

g2 <- ggplot(
  data = power |>
    dplyr::filter(id == "normal",  m == 2, delta < 5
                  # (delta < 5 & m == 2) |
                  # (delta < 13 & m == 5) |
                  #   m == 10
    ) |>
    dplyr::mutate(alternative = alt,
                  cens = factor(
                    dplyr::case_when((icens == "FALSE" & lcens == "0.0") ~ "none",
                                     (icens == "FALSE" & lcens == "0.2") ~ "left-censored",
                                     (icens == "TRUE" & lcens == "0.0") ~ "rounded",
                                     (icens == "TRUE" & lcens == "0.2") ~ "censored"))
    ),
  mapping = aes(
    x = delta,
    y = 100 * power,
    linetype = nobs,
    col = cens
  )
) +
  geom_hline(linetype = "dashed", alpha = 0.5, yintercept = 5) +
  geom_line() +
  labs(
    x = expression(delta),
    linetype = "sample size",
    color = "censoring",
    y = "power (%)",
    subtitle = "Generalized extreme value approximation to normal"
  ) +
  facet_wrap(
    facets = vars(alternative),
    labeller = label_both,
    scales = "free_x") +
  MetBrewer::scale_color_met_d(name = "Hiroshige") +
  theme_minimal() +
  theme(legend.position = 'bottom')

(g1 / g2) +
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')
ggsave(
  filename = "fig/power-curve-GEV-5.pdf",
  width = 10,
  height = 8,
  units = "in"
)
