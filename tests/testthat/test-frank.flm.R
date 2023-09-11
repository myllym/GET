test_that("frank.flm-methods", {
  nargs <- 3
  nobs <- 5
  nperm <- 22
  iset <- create_curve_set(list(obs=matrix(runif(nargs*nobs), nargs)))
  xset <- create_curve_set(list(obs=matrix(runif(nargs*nobs), nargs)))
  x <- data.frame(X=runif(nobs))

  # X = curve set
  f <- function(method) {
    set.seed(1)
    r <- frank.flm(nperm, Y ~ X, Y ~ 1, curve_sets = list(Y=iset, X=xset), typeone="fwer", savefuns=TRUE, method=method)
    attr(r, "runtime") <- NULL
    r
  }
  r1 <- f("lm")
  expect_equal(r1, f("ne"))
  expect_error(f("mlm"), "Curvesets in factors not allowed")
  expect_equal(r1, f("best"))

  # X = constant
  f <- function(method) {
    set.seed(1)
    r <- frank.flm(nperm, Y ~ X, Y ~ 1, curve_sets = list(Y=iset), factors=x, typeone="fwer", savefuns=TRUE, method=method)
    attr(r, "runtime") <- NULL
    r
  }
  r1 <- f("lm")
  expect_equal(r1, f("ne"))
  expect_equal(r1, f("best"))
  expect_equal(r1, f("mlm"))

  # X = curve_set + constant
  f <- function(method) {
    set.seed(1)
    r <- frank.flm(nperm, Y ~ X + Z, Y ~ 1, curve_sets = list(Y=iset, Z=xset), factors=x, typeone="fwer", savefuns=TRUE, method=method)
    attr(r, "runtime") <- NULL
    r
  }
  r1 <- f("lm")
  expect_equal(r1, f("ne"))
  expect_error(f("mlm"), "Curvesets in factors not allowed")
  expect_equal(r1, f("best"))
})

test_that("graph.flm-methods", {
  nargs <- 3
  nobs <- 5
  nperm <- 22
  iset <- create_curve_set(list(obs=matrix(runif(nargs*nobs), nargs)))
  xset <- create_curve_set(list(obs=matrix(runif(nargs*nobs), nargs)))
  x <- data.frame(X=runif(nobs), Z=runif(nobs))

  f <- function(method) {
    set.seed(1)
    r <- graph.flm(nperm, Y ~ X, Y ~ 1, curve_sets = list(Y=iset, X=xset), typeone="fwer", savefuns=TRUE, fast=method)
    attr(r, "runtime") <- NULL
    r
  }
  r1 <- f(FALSE)
  expect_equal(r1, f(TRUE))
  expect_equal(r1, f("mlm"))
  expect_equal(r1, f("ne"))
  expect_equal(r1, f("qr"))

  f <- function(method) {
    set.seed(1)
    r <- graph.flm(nperm, Y ~ X, Y ~ 1, curve_sets = list(Y=iset), factors=x, typeone="fwer", savefuns=TRUE, fast=method)
    attr(r, "runtime") <- NULL
    r
  }
  r1 <- f(FALSE)
  expect_equal(r1, f(TRUE))
  expect_equal(r1, f("mlm"))
  expect_equal(r1, f("ne"))
  expect_equal(r1, f("qr"))

  f <- function(method) {
    set.seed(1)
    r <- graph.flm(nperm, Y ~ X + Z, Y ~ 1, curve_sets = list(Y=iset), factors=x, typeone="fwer", savefuns=TRUE, fast=method)
    attr(r, "runtime") <- NULL
    r
  }
  r1 <- f(FALSE)
  expect_equal(r1, f(TRUE))
  expect_equal(r1, f("mlm"))
  expect_equal(r1, f("ne"))
  expect_equal(r1, f("qr"))

  f <- function(method) {
    set.seed(1)
    r <- graph.flm(nperm, Y ~ X + Z, Y ~ 1, curve_sets = list(Y=iset, Z=xset), factors=x, typeone="fwer", savefuns=TRUE, fast=method)
    attr(r, "runtime") <- NULL
    r
  }
  r1 <- f(FALSE)
  expect_equal(r1, f(TRUE))
  expect_equal(r1, f("mlm"))
  expect_equal(r1, f("ne"))
  expect_equal(r1, f("qr"))
})
