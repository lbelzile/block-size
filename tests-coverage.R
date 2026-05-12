# Check the coverage of binomial confidence intervals
setwd(this.path::here())
library(mev)
library(simsalapar)
source("helpers.R")
source("qqplot_NEW.R")


date <- 20250506
seed_init <- floor(runif(1)[1] * date)
get_alphastar_workhorse <- function(
  n,
  shape,
  m,
  rounding,
  leftcens,
  marginal,
  level = c(0.5, 0.8, 0.9, 0.95),
  nboot = 10000L,
  type = c("max",  "all")
) {
  xdat <- mev::build.blocks(
    mev::rgev(n = n, scale = 10, shape = shape),
    m = m
  )
  if (leftcens > 0) {
    lb <- mev::qgev(leftcens, scale = 10, shape = shape)
    type <- type[type != "range"]
  } else {
    lb <- NULL
  }
  if (rounding > 0) {
    xdat <- apply(xdat, 1:2, function(x) {
      round(x / rounding, 0) * rounding
    })
  }
  qqplot <- qqplot.blocksize(
    xdat = xdat,
    marginal = marginal,
    rounding = rounding,
    lb = lb,
    plot = FALSE,
    B = nboot,
    type = type,
    level = level[1],
    simult = TRUE
  )
  out <- apply(
    qqplot$gamma,
    2,
    quantile,
    probs = 1 - level,
    na.rm = TRUE
  )
  if (is.matrix(out)) {
    colnames(out) <- type
  } else if (is.vector(out)) {
    names(out) <- type
  }

  return(out)
}

get_coverage <- function(
  n,
  shape,
  m,
  rounding,
  leftcens,
  marginal,
  level,
  B = 2500L,
  type = c("max", "all")
) {
  if (leftcens > 0) {
    lb <- mev::qgev(leftcens, scale = 10, shape = shape)
  } else {
    lb <- NULL
  }
  if (is.vector(level)) {
    level <- matrix(level, nrow = 1)
  }
  cover <- total <- matrix(0, nrow = nrow(level), ncol = ncol(level))
  b <- 0
  while(b < B){
    xdat <- mev::build.blocks(
      mev::rgev(n = n, scale = 10, shape = shape),
      m = m
    )
    if (rounding > 0) {
      xdat <- apply(xdat, 1:2, function(x) {
        round(x / rounding, 0) * rounding
      })
    }
    for (i in seq_len(nrow(level))) {
      pp_boot <- try(
        qqplot.blocksize(
          xdat = xdat,
          simult = level[i, ],
          type = type,
          marginal = marginal,
          rounding = rounding,
          lb = lb,
          plot = FALSE
        ),
        silent = TRUE
      )
      if (!inherits(pp_boot, "try-error")) {
        b <- b + 1L
        for (j in seq_along(type)) {
          ypos <- pp_boot$plots[[j]]$y
          lower <- pp_boot$plots[[j]]$confint[, 1]
          upper <- pp_boot$plots[[j]]$confint[, 2]
          cover[i, j] <- cover[i, j] +
            isTRUE(all(ypos >= lower & ypos <= upper))
        }
      }
    }
  }
return(cover)
}

simu_fn_coverage <- function(sc){
  stopifnot(length(sc)==1L)
  shape_seq <- c(-0.25, 0, 0.25)
  n_seq <- c(25, 50, 100)
  m_seq <- c(2, 5)
  rounding_seq <- c(0, 1)
  leftcens_seq <- c(0, 0.2)
  marginal_seq <- c(FALSE, TRUE)

  cases <- expand.grid(
    shape = shape_seq,
    n = n_seq,
    m = m_seq,
    rounding = rounding_seq,
    leftcens = leftcens_seq,
    marginal = marginal_seq
  )
  alphastar <- get_alphastar_workhorse(n = cases$n[sc],
                          shape = cases$shape[sc],
                          m = cases$m[sc],
                          rounding = cases$rounding[sc],
                          leftcens = cases$leftcens[sc],
                          marginal = cases$marginal[sc],
                          level = c(0.5, 0.8, 0.9, 0.95),
                          type = c("max",  "all"),
                          nboot = 10000L)
  new <- as.list(cases[sc, ])
  new$level <- 1 - alphastar
  coverage <- do.call(
    get_coverage,
    args = new
  )
  list(alpha = alphastar, coverage = coverage)
}

nsc <- 144
nsim <- 5L
#Set list of variables for the simulation study
varList <- simsalapar::varlist(
  n.sim = list(type = "N", expr = quote(N), value = nsim),
  sc = list(type = "grid", value = 1:nsc)
)

out <- doClusterApply(
  vList = varList,
  doAL = FALSE,
  sfile = paste0("tests-coverage.rds"),
  cluster = parallel::makeCluster(40, type = "PSOCK"),
  doOne = simu_fn_coverage,
  check = FALSE,
  exports = ls()
)

if(FALSE){

shape_seq <- c(-0.25, 0, 0.25)
n_seq <- c(25, 50, 100)
m_seq <- c(2, 5)
rounding_seq <- c(0, 1)
leftcens_seq <- c(0, 0.2)
marginal_seq <- c(FALSE, TRUE)

cases <- expand.grid(
  shape = shape_seq,
  n = n_seq,
  m = m_seq,
  rounding = rounding_seq,
  leftcens = leftcens_seq,
  marginal = marginal_seq
)

res <- simsalapar::getArray(
  simsalapar::mkAL(
    x = out,
    vList = varList,
    repFirst = TRUE,
    check = FALSE),
  comp = "value")
fres <- array(0, dim = c(nsc, dim(res[1,1,1][[1]])))
for(j in seq_len(nsc)){
  for(i in seq_len(nsim)){
    fres[j,,] <- fres[j,,] + res[2,j,i][[1]]
  }
}
}
