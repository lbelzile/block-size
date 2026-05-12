setwd(this.path::here())

remotes::install_github("lbelzile/mev")
library(mev)
library(ggplot2)
library(patchwork)

shape_seq <- round(seq(-0.4, 0.4, by = 0.1), 2)
m_seq <- c(1:12)
info0 <- info1 <- array(
  dim = c(5, length(m_seq), length(shape_seq)),
  dimnames = list(
    "qty" = c("loc", "scale", "shape", "overall", "quant95"),
    m = m_seq,
    shape = shape_seq
  )
)

m.xi <- function(m, xi) {
  out <- log(m)
  if (abs(xi) > 1e-10){
    out <- (m^xi - 1) / xi
  }
  out
}

m1.xi <- function(m, xi) {
  out <- (log(m))^2 / 2
  if (abs(xi) > 1e-10){
    out <- (1 + xi * m^xi * log(m) - m^xi) / xi^2
  }
  out
}

Jacobian <- function(m, sigma, xi) {
  # Jacobian to transform information matrix for X_{(m)} from theta_m to theta
  a <- m.xi(m, xi)
  a1 <- m1.xi(m, xi)
  b <- m^xi
  b1 <- b * log(m)
  A <- matrix(0, 3, 3)
  A[1, ] <- c(1, 0, 0)
  A[2, ] <- c(a, b, 0)
  A[3, ] <- c(sigma * a1, sigma * b1, 1)
  A
}

ypfn <- function(par, p = 0.95) {
  mu <- par[1]
  sigma <- par[2]
  xi <- par[3]
  m <- 1 / (-log(p))
  c(1, m.xi(m, xi), sigma * m1.xi(m, xi))
}

mu <- 0
sigma <- 1 # results should not depend on location and scale, so include them to check for coding errors
p <- 0.95
for (i in seq_along(m_seq)) {
  for (j in seq_along(shape_seq)) {
    par1 <- c(mu, sigma, shape_seq[j])
    Fisherinfo1 <- mev::gev.infomat(par = par1, method = "exp", nobs = m_seq[i])
    par2 <- c(
      mu + sigma * m.xi(m = m_seq[i], xi = shape_seq[j]),
      sigma * m_seq[i]^shape_seq[j],
      shape_seq[j]
    )
    Fisherinfo2 <- mev::gev.infomat(par = par2, method = "exp", nobs = 1)
    A <- Jacobian(m = m_seq[i], sigma = sigma, xi = shape_seq[j])
    V <- solve(Fisherinfo1)
    Vm <- solve(A %*% Fisherinfo2 %*% t(A))
    J <- ypfn(par = c(mu, sigma, shape_seq[j]), p = p)
    Vp <- t(J) %*% V %*% J
    Vpm <- t(J) %*% Vm %*% J
    info0[, i, j] <- c(
      sqrt(diag(V) / diag(Vm)),
      (det(Vm) / det(V))^(1 / 3) / m_seq[i]^(2 * shape_seq[j] / 3 + 1),
      sqrt(Vp / Vpm)
    )
  }
}


# calculations for other return levels

yp_seq <- c(10, 20, 50, 100, 200)
rownames(info1) <- paste0("r[", yp_seq, "]")
p_seq <- 1 - 1 / yp_seq

for (i in seq_along(m_seq)) {
  for (j in seq_along(shape_seq)) {
    par1 <- c(mu, sigma, shape_seq[j])
    Fisherinfo1 <- mev::gev.infomat(par = par1, method = "exp", nobs = m_seq[i])
    par2 <- mev::maxstable(par1, m = m_seq[i])
    Fisherinfo2 <- mev::gev.infomat(par = par2, method = "exp", nobs = 1)
    A <- Jacobian(m = m_seq[i], sigma = sigma, xi = shape_seq[j])
    V <- solve(Fisherinfo1)
    Vm <- solve(A %*% Fisherinfo2 %*% t(A))
    for (k in seq_along(p_seq)) {
      J <- ypfn(par = c(mu, sigma, shape_seq[j]), p = p_seq[k])
      Vp <- t(J) %*% V %*% J
      Vpm <- t(J) %*% Vm %*% J
      info1[k, i, j] <- sqrt(Vp / Vpm)
    }
  }
}

qty_labs <- c(
  "loc" = "mu",
  "scale" = "sigma",
  "shape" = "xi",
  #"overall" = "overall",
  "quant95" = "r[20]")
plot_data <- array2DF(info0, responseName = "ratio") |>
  dplyr::filter(qty != "overall") |>
  dplyr::mutate(
    m = as.integer(m),
    qty = factor(qty, levels = names(qty_labs), labels = qty_labs))


ggplot(
  data = plot_data,
  mapping = aes(x = m, y = ratio, group = factor(shape), colour = factor(shape))) +
  geom_point() +
  geom_line() +
  facet_wrap(~qty, labeller = label_parsed, ) +
  scale_x_continuous(breaks = 2 * c(1:6), limits = c(1, 12)) +
  scale_y_continuous(breaks = 0.2 * c(0:5), limits = c(0, 1)) +
  MetBrewer::scale_color_met_d("Hiroshige") +
  labs(color = expression(xi),
       y = "") +
  theme_minimal()


ggsave(
  filename = "fig/Efficiencies.pdf",
  width = 10,
  height = 6,
  unit = 'in'
)

# quantile efficiences ==================================

plot_data_qt <- array2DF(info1[-2,,], responseName = "ratio") |>
    dplyr::mutate(
    m = as.integer(m),
    qty = factor(qty, ordered = TRUE, levels = rownames(info1))
)

ggplot(
  data = plot_data_qt,
  mapping = aes(x = m, y = ratio, group = factor(shape), colour = factor(shape))) +
  geom_point() +
  geom_line() +
  facet_wrap(~qty, labeller = label_parsed, ) +
  scale_x_continuous(breaks = 2 * c(1:6), limits = c(1, 12)) +
  scale_y_continuous(breaks = 0.2 * c(0:5), limits = c(0, 1)) +
  MetBrewer::scale_color_met_d("Hiroshige") +
  labs(color = expression(xi),
       y = "") +
  theme_minimal()


ggsave(
  filename = "fig/QuantileEfficiencies.pdf",
  width = 10,
  height = 6,
  unit = 'in'
)
