#' Fits an GPD to data
#'
#' @inheritParams FitGevFlex
#' @return The fitted parameters
#' @export
FitGpdFlex <-function (data, xpar, fpar, start, likelihood = "standard", ...) {
  
  flog.debug("Running FitGpdFlex")
  flog.debug("Version={%s}", paste0(packageVersion("evdflex")))
  
  if (likelihood == "standard") {
    flog.debug("Entering standard routine")
    nllGpd <- function(par) {
      pmat <- fpar(par, xpar)
      scale <- pmat$scale
      shape <- pmat$shape
      if (any(scale <= 0)) return(1e+20)
      exponential <- (abs(shape) < 1e-06)
      y <- data/scale
      z <- 1 + shape * y
      if (any(z <= 0, na.rm = TRUE)) return(1e+20)
      nll <- (shape + 1)/shape * log(z)
      nll[exponential] <- y[exponential]
      sum(nll + log(scale), na.rm = TRUE)
    }
  }
  call <- match.call()
  opt <- optim(start, nllGpd, ...)
  gpd <- fpar(opt$par, xpar)
  out <- list(estimate = opt$par, std.err = rep(NA, length(opt$par))
              ,cov = NULL, deviance = 2 * opt$value
              ,convergence = opt$convergence
              ,counts  = opt$counts
              ,message = opt$message
              ,loc     = gpd$loc
              ,scale   = gpd$scale
              ,shape   = gpd$shape)
  cmat <- try(solve(opt$hessian), TRUE)
  if (!inherits(cmat, "try-error")) {
    out$std.err <- sqrt(diag(cmat))
    out$cov <- cmat
  }
  
  structure(c(out, call = call), class = "gpd")
}
