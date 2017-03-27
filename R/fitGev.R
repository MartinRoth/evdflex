#' Fit flexible GEV
#'
#' @param data the data (in case of likelihood = "range" this should be a
#'   matrix or data.frame with two columns)
#' @param start the start values
#' @param fpar the parameter model
#' @param xpar the covariates used in fpar
#' @param std.err should standard errors be returned
#' @param ... additional parameters passed to the optimization
#' @param likelihood either standard or grouped
#' @export
FitGevFlex <- function(data, start, fpar, xpar, likelihood = "standard",
                       std.err = TRUE, ...) {

  flog.debug("Running FitGevFlex")
  flog.debug("Version={%s}", paste0(packageVersion("evdflex")))

  if (likelihood == "standard") {
    flog.debug("Entering standard routine")
    nll.gev <- function(par) {
      pmat  <- fpar(par, xpar)
      loc   <- pmat$loc
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
  } else if (likelihood == "grouped") {
    flog.debug("Entering grouped routine")
    countData <- nrow(data)
    nll.gev <- function(par) {

      pmat  <- AssertCorrectParameters(par, fpar, xpar, countData)
      loc   <- pmat$loc
      scale <- pmat$scale
      shape <- pmat$shape

      if(any(scale <= 0)) return(1e+20)

      gumbel <- (abs(shape) < 1e-06)

      y1 <- (data[ , 1] - loc) / scale
      y2 <- (data[ , 2] - loc) / scale

      nll <- numeric(countData)

      #nll for data2 == data1:
      isSame <- y1 == y2

      z1 <- 1 + shape * y1
      z2 <- 1 + shape * y2
      if(any(z1 <= 0, na.rm = TRUE)) return(1e+20)
      if(any(z2 <= 0, na.rm = TRUE)) return(1e+20)

      nll[isSame]          <- log(scale[isSame]) +
        (1 + 1 / shape[isSame]) * log(z1[isSame]) +
        z1[isSame]^(-1 / shape[isSame])
      nll[isSame & gumbel] <- log(scale[isSame & gumbel]) +
        y1[isSame & gumbel] + exp(-y1[isSame & gumbel])

      #nll for data2 <> data1:
      F1 <- exp(-(z1[!isSame])^(-1/shape[!isSame]))
      F2 <- exp(-(z2[!isSame])^(-1/shape[!isSame]))
      F1[gumbel[!isSame]] <- exp(-exp(-y1[!isSame]))
      F2[gumbel[!isSame]] <- exp(-exp(-y2[!isSame]))
      nll[!isSame]          <- -log(F2 - F1)
      nll[!isSame & gumbel] <- -log(F2[gumbel[!isSame]] - F1[gumbel[!isSame]])

      # nll <- c(nlleq, nllne)
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
  flog.debug("FitGevFlex run successfully!")
  structure(c(out, call = call), class = "evd")
}

AssertCorrectParameters <- function(par, fpar, xpar, N) {
  pmat  <- fpar(par, xpar)
  if (length(pmat$loc)   == 1) pmat$loc   <- rep(pmat$loc, N)
  else if (length(pmat$loc) != N)   stop("Loc should have length 1 or N")
  if (length(pmat$scale) == 1) pmat$scale <- rep(pmat$scale, N)
  else if (length(pmat$scale) != N) stop("Scale should have length 1 or N")
  if (length(pmat$shape) == 1) pmat$shape <- rep(pmat$shape, N)
  else if (length(pmat$shape) != N) stop("Shape should have length 1 or N")
  pmat
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


