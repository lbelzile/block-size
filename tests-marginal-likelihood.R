setwd(this.path::here())
library(mev)
library(ggplot2)

source("helpers.R")


m_seq <- c(2, 5, 10)
n_seq <- c(25L, 50L, 100L)
B <- 1000L
shape_seq <- seq(from = -0.3, to = 0.3, by = 0.05)
n_shape <- length(shape_seq)
coefs <- array(
  NA,
  dim = c(B, length(n_seq), length(shape_seq), length(m_seq), 3, 3)
)
dimnames(coefs) <- list(
  rep = 1:B,
  n = n_seq,
  shape = round(shape_seq, 2),
  m = m_seq,
  par = c("loc", "scale", "shape"),
  type = c("margc", "marg", "full")
)
for (b in seq_len(B)) {
  for (ns in seq_along(n_seq)) {
    for (s in seq_along(shape_seq)) {
      for (m in seq_along(m_seq)) {
        set.seed(b)
        xdat <- build.blocks(
          rgev(n = n_seq[ns] * m_seq[m], shape = shape_seq[s]),
          block = 1,
          m = m_seq[m]
        )
        coef_marg <- try(fit.gevblock(xdat, marginal = TRUE))
        if (!inherits(coef_marg, "try-error")) {
          coefs[b, ns, s, m, , 1] <- coef_marg
        }
        coef_marg2 <- try(fit.gevblock(
          xdat,
          marginal = TRUE,
          constraint = FALSE
        ))
        start <- NULL
        if (!inherits(coef_marg, "try-error")) {
          coefs[b, ns, s, m, , 2] <- coef_marg2
          start <- coef_marg2
        }
        coef_full <- try(fit.gevblock(xdat, marginal = FALSE, start = start))
        if (!inherits(coef_marg, "try-error")) {
          coefs[b, ns, s, m, , 3] <- coef_full
        }
      }
    }
  }
  # if (b %% 10 == 0) {
  #   print(b)
  #   save(coefs, file = "results/marginal_lik_test.RData")
  # }
}

save(coefs, file = ".results/marginal_lik_test.RData")

coef_df <- array2DF(coefs, responseName = "coef") |>
  dplyr::mutate(
    n = factor(as.integer(n)),
    shape = as.numeric(shape),
    m = factor(as.integer(m)),
    par = factor(par),
    type = factor(type, labels = c("full", "marginal", "marginal + constraint"))
  )

ggplot(
  data = coef_df |>
    dplyr::filter(
      par == "shape",
      shape %in% c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),
      n %in% c("25", "100"),
      m %in% c("2", "10")
    ),
  mapping = aes(
    x = factor(shape),
    y = coef - shape,
    col = type,
    fill = type,
    group = type
  )
) +
  ggdist::stat_pointinterval(position = position_dodge()) +
  facet_wrap(
    vars(n, m),
    labeller = labeller(
      n = ~ paste0("n=", .),
      m = ~ paste0("m=", .),
      .multi_line = FALSE
    ),
    scales = "free_y"
  ) +
  labs(x = "shape", y = "bias") +
  MetBrewer::scale_color_met_d("Hiroshige") +
  MetBrewer::scale_fill_met_d("Hiroshige") +
  theme_classic() +
  theme(legend.position = "bottom")
