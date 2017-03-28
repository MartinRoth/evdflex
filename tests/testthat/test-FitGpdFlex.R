context("FitGpdFlex")

library(futile.logger)
flog.threshold(DEBUG)
flog.appender(appender.file('gpdflex.log'))

library(data.table)

testData <- fread("testData.txt")

testData <- as.data.frame(testData[, c(2,3), with = FALSE])

testData <- subset(testData, V2 >= 10)
testData <- testData - 10

test_that("regression test (gpd) - normal fit", {

  start <- c(2.3, 0)

  xpar <- list(N = nrow(testData))
  
  fpar <- function(p, xpar) {
    loc   <- 10
    scale <- p[1]
    shape <- p[2]
    list(loc = loc, scale = scale, shape = shape)
  }

  tmp <- FitGpdFlex(data = testData$V2,
                    xpar = xpar,
                    fpar = fpar,
                    start = start
                    )

  expect_equal_to_reference(tmp, file = "./referenceOutput/normalFitGpd.rds")
  
  start <- c(2.3)
  
  xpar <- list(N = nrow(testData))
  
  fpar <- function(p, xpar) {
    loc   <- 10
    scale <- p[1]
    shape <- 0
    list(loc = loc, scale = scale, shape = shape)
  }
  
  tmp <- FitGpdFlex(data = testData$V2,
                    xpar = xpar,
                    fpar = fpar,
                    start = start,
                    method = "Brent",
                    lower = 0.1,
                    upper = 15)
  
  expect_equal_to_reference(tmp, file = "./referenceOutput/normalFitExponential.rds")
  
})

test_that("duplicates only", {
  start <- c(2.3, 0)
  
  xpar <- list(N = nrow(testData))
  
  fpar <- function(p, xpar) {
    loc   <- 10
    scale <- p[1]
    shape <- p[2]
    list(loc = loc, scale = scale, shape = shape)
  }
  
  tmp1 <- FitGpdFlex(data = testData$V2,
                     xpar = xpar,
                     fpar = fpar,
                     start = start
  )
  tmp1$call <- NULL
  
  tmp2 <- FitGpdFlex(data = data.frame(testData$V2, testData$V2),
                    xpar = xpar,
                    fpar = fpar,
                    start = start,
                    likelihood = "grouped")
  tmp2$call <- NULL
  
  expect_equal(tmp1, tmp2)
  
})

test_that("regression test - grouped fit", {

  start <- c(2.3, 0)
  
  xpar <- list(N = nrow(testData))
  
  fpar <- function(p, xpar) {
    loc   <- 10
    scale <- p[1]
    shape <- p[2]
    list(loc = loc, scale = scale, shape = shape)
  }
  
  tmp <- FitGpdFlex(data = testData,
                    xpar = xpar,
                    fpar = fpar,
                    start = start,
                    likelihood = "grouped")
  
  start <- c(2.3)
  
  xpar <- list(N = 1)
  
  fpar <- function(p, xpar) {
    loc   <- 10
    scale <- p[1]
    shape <- 0
    list(loc = loc, scale = scale, shape = shape)
  }
  
  tmp <- FitGpdFlex(data = data.frame(2, 4),
                    xpar = xpar,
                    fpar = fpar,
                    start = start,
                    likelihood = "grouped",
                    method = "Brent", lower = 0.1, upper = 10)
  
  expect_equal(tmp$scale, 2.88539, tolerance = 1e-6)
  #expect_equal_to_reference(tmp, file = "./referenceOutput/tinyFitGpd.rds")

#   tmp <- FitGevFlex(testData[, c(2, 3), with = FALSE], start = start,
#                     fpar = fpar, likelihood = "grouped")
# 
#   expect_equal_to_reference(tmp, file = "./referenceOutput/groupedFit.rds")
# 
#   xpar <- list(N = nrow(testData))
# 
#   fpar2 <- function(p, xpar) {
#     loc   <- rep(p[1], xpar$N)
#     scale <- rep(p[2], xpar$N)
#     shape <- rep(p[3], xpar$N)
#     list(loc = loc, scale = scale, shape = shape)
#   }
# 
#   tmp2 <- FitGevFlex(testData[, c(2, 3), with = FALSE], start = start,
#                      fpar = fpar2, xpar=xpar,
#                      likelihood = "grouped")
#   expect_equal(tmp2$estimate, tmp$estimate)
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
