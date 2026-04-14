setwd(this.path::here())

remotes::install_github("paulnorthrop/evils")
library(ggplot2)
library(mev)
library(evils)
library(numDeriv)

# Cheap estimation  of the Fisher information matrix via Monte Carlo

# Values for the shape
shape_seq <- round(seq(-0.3, 0.5, by = 0.1), 10)
nxi <- length(shape_seq)
m_seq <- c(2L, 5L, 10L) # block size
nm <- length(m_seq)


# Marginal log-likelihood for m-1 lowest order statistics from block of m
ll_marg_fn <- function(par, xdat) {
  m <- ncol(xdat)
  sum(evils::dGEV(
    x = c(xdat[, -m]),
    loc = par[1],
    scale = par[2],
    shape = par[3],
    log = TRUE
  )) +
    sum(evils::pGEV(
      q = c(xdat[, m - 1]),
      loc = par[1],
      scale = par[2],
      shape = par[3],
      lower.tail = FALSE,
      log.p = TRUE
    ))
}


n <- 1000L # sample size
B <- 1000L # number of Monte Carlo repetitions


info_marg <- array(dim = c(length(m_seq), length(shape_seq), 3, 3, B))
for (i in seq_along(m_seq)) {
  for (j in seq_along(shape_seq)) {
    for (b in 1:B) {
      set.seed(b)
      xdat <- matrix(
        mev::rgev(
          n = n * m_seq[i],
          loc = 0,
          scale = 1,
          shape = shape_seq[j]
        ),
        ncol = m_seq[i],
        nrow = n
      )
      xdat <- t(apply(xdat, 1, sort))
      # Compute information as square of score
      info_marg[i, j, , , b] <- tcrossprod(numDeriv::grad(
        func = ll_marg_fn,
        x = c(0, 1, shape_seq[j]),
        xdat = xdat
      )) /
        n
    }
  }
}

# Calculate average
info_marg_mean <- apply(info_marg, 1:4, mean, na.rm = TRUE)
# Check that precision is good enough for reporting
info_marg_sd <- apply(info_marg, 1:4, sd, na.rm = TRUE) / sqrt(B)

# Extract efficiency - cubic root of ratio of determinants
eff <- matrix(nrow = nm, ncol = nxi)
for (i in seq_along(m_seq)) {
  for (j in seq_along(shape_seq)) {
    eff[i, j] <- exp(
      (determinant(info_marg_mean[i, j, , ])$modulus -
        determinant(mev::gev.infomat(
          par = c(0, 1, shape_seq[j]),
          nobs = m_seq[i],
          method = "exp"
        ))$modulus) /
        3
    )
  }
}
dimnames(eff) <- list(m = m_seq, shape = round(shape_seq, 2))

# Plot efficiency as a function of the shape
matplot(
  shape_seq,
  t(eff),
  ylab = "efficiency",
  xlab = expression(xi),
  bty = "l",
  type = "b"
)

# Make into a data frame and use grammar of graphics
eff_df <- array2DF(x = eff, responseName = "efficiency") |>
  dplyr::mutate(
    shape = as.numeric(shape),
    m = factor(m, levels = c(2, 5, 10), labels = c("2", "5", "10"))
  )

ggplot(
  data = eff_df,
  mapping = aes(x = shape, col = m, linetype = m, y = efficiency)
) +
  geom_point() +
  geom_line() +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(),
    labels = c("0", "0.25", "0.5", "0.75", "1")
  ) +
  MetBrewer::scale_color_met_d("Hiroshige") +
  theme_classic() +
  theme(legend.position = "inside", legend.position.inside = c(0.9, 0.1))
# Save results for future use
save(eff_df, file = "results/efficiency.RData")
