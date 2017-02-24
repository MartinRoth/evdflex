context("FitGevFlex")

## TODO: Rename context
## TODO: Add more tests

library(futile.logger)
flog.threshold(DEBUG)
flog.appender(appender.file('gevflex.log'))

library(data.table)

testData <- fread("./testData.txt")

test_that("regression test - normal fit", {

  start <- c(7.5, 2.3, 0)

  fpar <- function(p, xpar) {
    loc   <- p[1]
    scale <- p[2]
    shape <- p[3]
    list(loc = loc, scale = scale, shape = shape)
  }

  tmp <- FitGevFlex(testData[, 2, with = FALSE], start = start, fpar = fpar)

  expect_equal_to_reference(tmp, file = "./referenceOutput/normalFit.rds")
})

test_that("regression test - range fit", {

  start <- c(7.5, 2.3, 0)

  fpar <- function(p, xpar) {
    loc   <- p[1]
    scale <- p[2]
    shape <- p[3]
    list(loc = loc, scale = scale, shape = shape)
  }

  tmp <- FitGevFlex(testData[, c(2, 3), with = FALSE], start = start,
                    fpar = fpar, likelihood = "range")

  expect_equal_to_reference(tmp, file = "./referenceOutput/rangeFit.rds")

  xpar <- list(N = nrow(testData))

  fpar2 <- function(p, xpar) {
    loc   <- rep(p[1], xpar$N)
    scale <- rep(p[2], xpar$N)
    shape <- rep(p[3], xpar$N)
    list(loc = loc, scale = scale, shape = shape)
  }

  tmp2 <- FitGevFlex(testData[, c(2, 3), with = FALSE], start = start,
                     fpar = fpar2, xpar=xpar,
                     likelihood = "range")
  expect_equal(tmp2$estimate, tmp$estimate)
})

# test_that("regression test - one dimensional optim", {
#
#   start <- 7.5
#
#   fpar <- function(p, xpar) {
#     loc   <- p[1]
#     scale <- 2.5286
#     shape <- 0.0934
#     list(loc = loc, scale = scale, shape = shape)
#   }
#
#   tmp <- fgev.flex_Brent(testData[, 2, with = FALSE], start = start,
#                          fpar = fpar)
#
#   expect_equal_to_reference(tmp, file = "./referenceOutput/oneDimFit.rds")
# })
