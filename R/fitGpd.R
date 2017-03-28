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
  } else if (likelihood == "grouped") {
    flog.debug("Entering grouped routine")
    countData <- nrow(data)
    nllGpd <- function(par) {
      
      pmat  <- AssertCorrectParameters(par, fpar, xpar, countData)
      loc   <- pmat$loc
      scale <- pmat$scale
      shape <- pmat$shape
      
      if(any(scale <= 0)) return(1e+20)
      
      exponential <- (abs(shape) < 1e-06)
      
      # y1 <- (data[ , 1] - loc) / scale
      # y2 <- (data[ , 2] - loc) / scale
      y1 <- data[, 1] / scale
      y2 <- data[, 2] / scale
      
      nll <- numeric(countData)
      
      #nll for data2 == data1:
      isSame <- y1 == y2
      
      z1 <- 1 + shape * y1
      z2 <- 1 + shape * y2
      if(any(z1 <= 0, na.rm = TRUE)) return(1e+20)
      if(any(z2 <= 0, na.rm = TRUE)) return(1e+20)
      
      nll[isSame] <- (shape[isSame] + 1) / shape[isSame] * log(z1[isSame])
      nll[isSame & exponential] <- y1[isSame & exponential]
      nll[isSame] <- nll[isSame] + log(scale[isSame])
      
      #nll for data2 <> data1:
      F1 <- -z1[!isSame]^(-1/shape[!isSame])
      F2 <- -z2[!isSame]^(-1/shape[!isSame])
      F1[exponential[!isSame]] <- -exp(-y1[!isSame]) 
      F2[exponential[!isSame]] <- -exp(-y2[!isSame]) 
      nll[!isSame]          <- -log(F2 - F1)
      nll[!isSame & exponential] <- -log(F2[exponential[!isSame]] - F1[exponential[!isSame]])
      
      sum(nll, na.rm = TRUE)
    }  
  } else {
    stop("Likelihood not determined.")
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
