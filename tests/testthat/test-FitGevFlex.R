context("FitGevFlex")

## TODO: Rename context
## TODO: Add more tests

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

  tmp <- fgev.flex_range(testData[, 2, with = FALSE],
                         testData[, 3, with = FALSE],
                         start = start, fpar = fpar)

  expect_equal_to_reference(tmp, file = "./referenceOutput/rangeFit.rds")
})
