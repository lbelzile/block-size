#' Simulation from alternative for the block size simulation study
#'
#' Given two sets of parameters for the generalized extreme value
#' distribution, generate \code{nm} observations from the model,
#' where the largest observation of \code{m} in each of the
#' \code{n} block is drawn from the model from \code{coef2},
#' which is left-truncated at the penultimate value
#' @param n sample size
#' @param m number of order statistics
#' @param coef1 vector of c GEV parameters for the (\code{m}-1) smallest order statistics
#' @param coef2 vector of GEV parameters for the largest of \code{m} observations
#' @param seed optional integer to set the seed (if not \code{NULL})
#' @return an \code{n} by \code{m} matrix of observations, ordered by rows
#' @examples
#' delta <- 0.1
#' coef1 <- c(0,1,0.1)
#' coef2 <- c(coef1[1] + delta, coef1[2] * exp(-delta / 10), coef1[3])
#' sim_alt(n = 100, m = 4, coef1 = coef1, coef2 = coef2)
sim_alt <- function(n, m, coef1, coef2, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(as.integer(seed[1]))
  }
  coef1 <- as.numeric(coef1)[1:3]
  coef2 <- as.numeric(coef2)[1:3]
  xdat <- matrix(
    mev::rgev(
      n = m * n,
      loc = coef1[1],
      scale = coef1[2],
      shape = coef1[3]
    ),
    nrow = n,
    ncol = m
  )
  xdat <- t(apply(xdat, 1, sort))
  # Coefficients for the block max
  if (max(abs(coef1 - coef2)) > 1e-8) {
    # The largest observation is truncated GEV
    Fa <- mev::pgev(
      q = xdat[, m - 1],
      loc = coef2[1],
      scale = coef2[2],
      shape = coef2[3]
    )
    xdat[, m] <- mev::qgev(
      p = Fa + runif(n) * (1 - Fa),
      loc = coef2[1],
      scale = coef2[2],
      shape = coef2[3]
    )
  }
  return(xdat)
}


#' Simulate block maxima from a distribution
#'
#' Simulate from the distribution of a \code{m} observation maximum
#' @param n number of observations
#' @param m block size
#' @param ltrunc lower truncation bound (if none, use \code{NULL}); only supoorted by \code{method="inversion"}
#' @param args list containing the additional parameters of the distribution, e.g., \code{shape1} and \code{shape2} for \code{dist="beta"}
#' @param dist string indicating the distribution
#' @param method string; either \code{inversion} of the CDF, or \code{max} for direct simulation.
#' @param vector of length \code{n}
#' @examples
#' rdistmax(n = 100, dist = "norm", args = list(sd = 10), m = 30)
#' rdistmax(n = 100, dist = "norm", args = list(sd = 10), ltrunc = 25, m = 30, method = "max")
rdistmax <- function(
  n,
  dist,
  m,
  ltrunc = NULL,
  args = NULL,
  method = c("inversion", "max")
) {
  method <- match.arg(method)
  if (is.null(args)) {
    args <- list()
  }
  if (!is.null(ltrunc)) {
    if (length(ltrunc) %in% c(1L, n)) {
      ltrunc <- rep(ltrunc, length.out = n)
    }
  }
  if (method == "inversion") {
    if (is.null(ltrunc)) {
      args$p <- -rexp(n) / m
    } else {
      args$q <- ltrunc
      args$log.p <- TRUE
      F_a <- m * do.call(paste0("p", dist), args = args)
      args$q <- NULL
      args$p <- log(exp(F_a) + (1 - exp(F_a)) * runif(n)) / m
    }
    args$log.p <- TRUE
    return(do.call(paste0("q", dist), args = args))
  } else {
    args$n <- m * n
    out <- apply(
      matrix(do.call(paste0("r", dist), args = args), nrow = n, ncol = m),
      1,
      max
    )
    if (is.null(ltrunc)) {
      return(out)
    }
    nfail <- sum(out < ltrunc)
    while (nfail > 0) {
      args$n <- nfail * m
      out[out < ltrunc] <- apply(
        matrix(do.call(paste0("r", dist), args = args), nrow = nfail, ncol = m),
        1,
        max
      )
      nfail <- sum(out < ltrunc)
    }
    return(out)
  }
}

#' Simulate from max domain of attraction
#'
#' Simulate \code{n} vectors of order statistics, with
#' each simulate from a distribution \code{dist},
#' each value corresponding to the maximum of \code{m0}
#' (or \code{m1}) independent observations.
#'
#' @param n number of observations
#' @param m block size
#' @param m0 block size of original distribution for null
#' @param m1 block size of original distribution of \code{m}th order statistic
#' @param args list containing the additional parameters of the distribution, e.g., \code{shape1} and \code{shape2} for \code{dist="beta"}
#' @param dist string indicating the distribution
#' @param seed optional integer to set the seed (if not \code{NULL})
#' @return an \code{n} by \code{m} matrix of observations, ordered by rows
sim_alt_mda <- function(n, m, m0, m1, dist, args = NULL, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  m <- as.integer(m)[1]
  n <- as.integer(n)[1]
  xdat <- matrix(
    data = rdistmax(n = n * m, m = m0, dist = dist, args = args),
    nrow = n,
    ncol = m
  )
  xdat <- t(apply(xdat, 1, sort))
  # Coefficients for the block max
  if (abs(m0 - m1) > 1e-8) {
    xdat[, m] <- rdistmax(
      n = n,
      m = m1,
      ltrunc = xdat[, m - 1],
      dist = dist,
      args = args
    )
  }
  return(xdat)
}


#' Power study from stationary distribution
#'
#' This function generate samples of observations, run
#' likelihood ratio tests for serially dependent data from
#' max first-order autoregressive process with tail index \code{theta}
#' @param nobs number of vectors of size m to generate
#' @param m block size
#' @param theta extremal index, between (0,1]
#' @param shape non-negative scalar shape parameter for the GEV
#' @param alt integer 1, 2 or 3, indicating the alternative distribution
#' @param block block size over which to compute maxima
#' @return a vector of p-values of the same length as \code{alt}
simu_fn_st <- function(
  nobs,
  m,
  theta = 1,
  shape,
  alt = 1:3,
  block = 1L
) {
  block <- as.integer(block)
  stopifnot(block >= 1L, shape >= 0)
  xdat <- array(
    mev::rmar1(n = nobs * m * block, theta = theta, shape = shape),
    dim = c(block, nobs, m)
  )
  xdat <- apply(xdat, 2:3, max)
  xdat <- t(apply(xdat, 1, sort))
  pvals <- try(test.blocksize(xdat, alt = alt), silent = TRUE)
  if (inherits(pvals, "try-error")) {
    pvals <- rep(NA, length(alt))
  }
  if (is.data.frame(pvals)) {
    pvals <- pvals$pval
  }
  return(pvals)
}


#' Power study from iid GEV data
#'
#' This function generate samples of observations, run
#' likelihood ratio tests for serially dependent data from
#' max first-order autoregressive process with tail index \code{theta}
#' @param nobs number of vectors of size m to generate
#' @param m block size
#' @param delta vector of positive (or larger than 1) indicating departure from max-stability
#' @param id integer indicator for scenario, either \code{1} for GEV(0,1,0.1) versus GEV(\eqn{\delta, \exp(-\delta/10), 0.1}), \code{2} for the GEV penultimate approximation to blocks of size 30 from a Weibull(1, 0.8) or \code{3} for GEV penultimate approximation to blocks of size 30 from standard normal.
#' @param alt integer 1, 2 or 3, indicating the alternative distribution
#' @return a vector of p-values of the same length as \code{alt}
simu_fn <- function(id, delta, nobs, m, alt) {
  stopifnot(isTRUE(all(alt %in% 1:4)))
  # Simulate max-stable data from GEV
  if (id == 1L) {
    pars1 <- c(0, 1, 0.1)
    pars2 <- c(delta, exp(-delta / 10), 0.1)
    xdat <- sim_alt(n = nobs, m = m, coef1 = pars1, coef2 = pars2)
  } else if (id == 2L) {
    # Simulate max-stable data from GEV
    pars1 <- mev::penultimate(
      family = 'weibull',
      shape = 0.8,
      m = 30,
      method = "bm",
      returnList = FALSE
    )[1, 1:3]
    pars2 <- mev::penultimate(
      family = 'weibull',
      shape = 0.8,
      m = 30 * delta,
      method = "bm",
      returnList = FALSE
    )[1, 1:3]
    xdat <- sim_alt(
      n = nobs,
      m = m,
      coef1 = pars1,
      coef2 = pars2
    )
  } else if (id == 3L) {
    pars1 <- mev::penultimate(
      family = 'norm',
      m = 30,
      method = "bm",
      returnList = FALSE
    )[1, 1:3]
    pars2 <- mev::penultimate(
      family = 'norm',
      m = 30 * delta,
      method = "bm",
      returnList = FALSE
    )[1, 1:3]
    xdat <- sim_alt(
      n = nobs,
      m = m,
      coef1 = pars1,
      coef2 = pars2
    )
  }

  pvals <- try(test.blocksize(xdat, alt = alt), silent = TRUE)
  if (inherits(pvals, "try-error")) {
    pvals <- rep(NA, length(alt))
  }
  if (is.data.frame(pvals)) {
    pvals <- pvals$pval
  }
  return(pvals)
}


#' Power study from data from the max domain of attraction of a GEV
#'
#' This function generate samples of observations, run
#' likelihood ratio tests for serially dependent data from
#' max first-order autoregressive process with tail index \code{theta}
#' @param nobs number of vectors of size m to generate
#' @param m block size
#' @param delta vector of positive (or larger than 1) indicating departure from max-stability
#' @param id integer indicator for scenario, either \code{2} for the GEV penultimate approximation to blocks of size \code{m0} from a Weibull(1, 0.8) or \code{3} for GEV penultimate approximation to blocks of size \code{m0} from standard normal.
#' @param alt integer 1, 2 or 3, indicating the alternative distribution
#' @param m0 size of block
#' @return a vector of p-values of the same length as \code{alt}
simu_fn_mda <- function(id, delta, nobs, m, alt, m0 = 30) {
  stopifnot(delta >= 1)
  stopifnot(isTRUE(all(alt %in% 1:3)), isTRUE(all(id %in% 2:3)))
  # Simulate max-stable data from GEV
  if (id == 2L) {
    dist = "weibull"
    args <- list(shape = 0.8)
  } else if (id == 3L) {
    dist <- "norm"
    args <- NULL
  }
  xdat <- sim_alt_mda(
    n = nobs,
    m = m,
    m0 = m0,
    m1 = m0 * delta,
    dist = dist,
    args = args
  )
  pvals <- try(get_pvals(xdat, alt = alt), silent = TRUE)
  if (is.data.frame(pvals)) {
    pvals <- pvals$pval
  }
  if (inherits(pvals, "try-error")) {
    pvals <- rep(NA, length(alt))
  }
  return(pvals)
}


# ggplot2 function for objects of class `mev_plot_blocksize`
# produced by the function `qqplot.blocksize` in package `mev`.
#' @export
autoplot.mev_plot_blocksize <- function(
  x,
  type = c("max", "spacing", "all"),
  ...
) {
  type <- match.arg(type)
  type_plot <- lapply(x$plots, function(x) {
    x$type
  })
  tind <- which(type_plot == type)

  df <- with(
    x$plots[[tind]],
    data.frame(
      x = x,
      y = y,
      jtlower = confint$simultaneous[, 1],
      jtupper = confint$simultaneous[, 2],
      ptlower = confint$pointwise[, 1],
      ptupper = confint$pointwise[, 2]
    )
  )
  gg1 <- ggplot(data = df) +
    geom_abline(intercept = 0, slope = 1) +
    geom_segment(
      alpha = 0.1,
      col = "grey90",
      mapping = aes(x = x, y = jtlower, yend = jtupper)
    ) +
    geom_segment(
      alpha = 0.2,
      col = "grey70",
      mapping = aes(x = x, y = ptlower, yend = ptupper)
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      expand = expansion(add = c(0, 0.01)),
      breaks = seq(0, 1, by = 0.25),
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(add = c(0, 0.01)),
      breaks = seq(0, 1, by = 0.25),
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) +
    geom_point(
      mapping = aes(
        x = x,
        y = y,
        col = factor(ifelse(y < jtlower | y > jtupper, "outside", "inside"))
      )
    ) +
    scale_color_manual(values = c("black", "#e76254")) +
    labs(x = "theoretical positions", y = "empirical positions") +
    theme_classic() +
    theme(legend.position = "none")
  return(gg1)
}
