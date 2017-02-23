#' Fit flexible GEV
#'
#' @param data the data (in case of likelihood = "range" this should be a
#'   matrix or data.frame with two columns)
#' @param start the start values
#' @param fpar the parameter model
#' @param xpar the covariates used in fpar
#' @param std.err should standard errors be returned
#' @param ... additional parameters passed to the optimization
#' @param likelihood either standard or range
#' @export
FitGevFlex <- function(data, start, fpar, xpar, likelihood = "standard",
                       std.err = TRUE, ...) {

  if (likelihood == "standard") {
    nll.gev <- function(par) {
      pmat <- fpar(par, xpar)
      loc <- pmat$loc
      scale <- pmat$scale
      shape <- pmat$shape
      if(any(scale <= 0)) return(1e+20)
      gumbel <- (abs(shape) < 1e-06)
      y <- (data - loc) / scale
      z <- 1 + shape * y
      if(any(z <= 0, na.rm = TRUE)) return(1e+20)
      nll <- (1 + 1 / shape) * log(z) + z^(-1 / shape)
      nll[gumbel] <- y[gumbel] + exp(-y[gumbel])
      sum(nll + log(scale), na.rm = TRUE)
    }
  } else if (likelihood == "range") {
    nll.gev <- function(par) {
      pmat <- fpar(par, xpar)
      loc <- pmat$loc
      scale <- pmat$scale
      shape <- pmat$shape
      if(any(scale <= 0)) return(1e+20)
      gumbel <- (abs(shape) < 1e-06)
      #   sum(nll + log(scale), na.rm = TRUE)

      y1 <- (data[ , 1] - loc) / scale
      y2 <- (data[ , 2] - loc) / scale

      #nll for data2 == data1:
      y1eq <- y1[y1 == y2]
      z1eq <- 1 + shape * y1eq
      if(any(z1eq <= 0, na.rm = TRUE)) return(1e+20)
      nlleq         <- log(scale) + (1 + 1 / shape) * log(z1eq) + z1eq^(-1 / shape)
      nlleq[gumbel] <- log(scale) + y1eq[gumbel] + exp(-y1eq[gumbel])

      #nll for data2 <> data1:
      y1ne <- y1[y1 != y2]
      y2ne <- y2[y1 != y2]
      z1ne <- 1 + shape * y1ne
      z2ne <- 1 + shape * y2ne
      if(any(z1ne <= 0, na.rm = TRUE)) return(1e+20)
      if(any(z2ne <= 0, na.rm = TRUE)) return(1e+20)
      F1 <- exp(-(z1ne)^(-1/shape))
      F2 <- exp(-(z2ne)^(-1/shape))
      F1[gumbel] <- exp(-exp(-y1ne))
      F2[gumbel] <- exp(-exp(-y2ne))
      nllne         <- -log(F2 - F1)
      nllne[gumbel] <- -log(F2[gumbel] - F1[gumbel])

      nll <- c(nlleq, nllne)
      sum(nll, na.rm = TRUE)
    }
  } else {
    stop("Likelihood not determined.")
  }
  call <- match.call()
  opt <- optim(start, nll.gev, hessian = std.err, ...)
  gev <- fpar(opt$par, xpar)
  out <- list(estimate = opt$par, std.err = rep(NA, length(opt$par)),
              cov = NULL, deviance = 2 * opt$value,
              convergence = opt$convergence, counts = opt$counts,
              message = opt$message,
              loc = gev$loc, scale = gev$scale, shape = gev$shape)
  cmat <- try(solve(opt$hessian), TRUE)
  if(!inherits(cmat, "try-error")) {
    out$std.err <- sqrt(diag(cmat))
    out$cov <- cmat
  }
  structure(c(out, call = call), class = "evd")
}




####################

#Added to solve one-dimensional fits
#Jules Beersma
#20150326

fgev.flex_Brent <- function(data, start, fpar, xpar, std.err = TRUE, ...) {
  nll.gev <- function(par) {
    pmat <- fpar(par, xpar)
    loc <- pmat$loc
    scale <- pmat$scale
    shape <- pmat$shape
    if(any(scale <= 0)) return(1e+20)
    gumbel <- (abs(shape) < 1e-06)
    y <- (data - loc) / scale
    z <- 1 + shape * y
    if(any(z <= 0, na.rm = TRUE)) return(1e+20)
    nll <- (1 + 1 / shape) * log(z) + z^(-1 / shape)
    nll[gumbel] <- y[gumbel] + exp(-y[gumbel])
    sum(nll + log(scale), na.rm = TRUE)
  }
  call <- match.call()
  #  opt <- optim(start, nll.gev, hessian = std.err, ...)
  opt <- optim(start, nll.gev, method = "Brent", lower = 5, upper = 200, hessian = std.err, ...)
  #  opt <- optimize(nll.gev, interval = c(5, 200))
  gev <- fpar(opt$par, xpar)
  out <- list(estimate = opt$par, std.err = rep(NA, length(opt$par)), cov = NULL, deviance = 2 * opt$value, convergence = opt$convergence, counts = opt$counts, message = opt$message, loc = gev$loc, scale = gev$scale, shape = gev$shape)
  cmat <- try(solve(opt$hessian), TRUE)
  if(!inherits(cmat, "try-error")) {
    out$std.err <- sqrt(diag(cmat))
    out$cov <- cmat
  }
  structure(c(out, call = call), class = "evd")
}


####################


