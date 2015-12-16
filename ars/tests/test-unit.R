library(testthat)
library(ars)
library(numDeriv)
source('KS_stats.R')
source('function_for_test.R')

test_that("Testing ConstructZ:", {
  print("Unit testing ConstructZ function:")
  expect_equal(ConstructZ(0, dnorm, -2,2),c(-2,2), tolerance = .001)
  expect_equal(ConstructZ(c(-2,2), dnorm,-2,2),c(-2,0,2))
  expect_equal (ConstructZ(c(-100,2,100), dnorm,-2,2),c(-2,2.5,2.5,2), tolerance = .001)
  print("passed ConstructZ unit test")
})

#the following test works when run locally, but doesn't allow for package to be built
#properly displays warning on non-log-concave functions as inputs
#test_that("Unit testing log-concavity function:",{
  #h <- function(x) x^2
  #print("Unit testing logconc function:")
  #expect_warning(logconc(h,-1,1))
  #print("Passed non-log concave test")
#})

test_that("Unit testing log-concavity function:",{
  h <- dnorm
  print("Unit testing logconc function:")
  print(expect_null(logconc(h,-4,4)))
  print("Passed log concave test")
})

test_that("Unit testing findmax function:",{
  #standard case with one max
  h <- dnorm
  lb <- -1
  ub <- 1
  print("Unit testing findmax function:")
  print(expect_equal(findmax(h,-1,1),c(-.002,0,.002),tolerance = .00001))
  print("Passed findmax test")

  #tricky case with 2 max points
  h <- function(x) x^2
  lb <- -1
  ub <- 1
  print("Unit testing findmax function:")
  #tolerance here is delta in function
  print(expect_equal(findmax(h,-1,1),c(-1,-1),tolerance = .01))
  print("Passed findmax test")
})

test_that("Testing utest function for uniform distribution:", {
  print("Unit testing utest function:")
  g <- dunif
  expect_equal(utest(g,0,1),TRUE)
  print("passed uniform case utest unit test")
})

test_that("Testing utest function for uniform distribution:", {
  print("Unit testing utest function:")
  g <- dnorm
  expect_equal(utest(g,-1,1),FALSE)
  print("passed non-uniform case utest unit test")
})

#tesetthat

#library(testthat)
#source('main.R')
#Test for upper
#Around 0, the tagent line has a slope of zero ,so the value near 0 should all be equal
g<- dnorm
test_that("Upper", {
  print("Testing Upper function now:")
  X_k <- c(-1,0,1)
  lb <- -1
  ub <- 1
  H_k <- h(X_k)
  dH_k <- grad(h,X_k)
  Z_k <- ConstructZ(X_k, h, lb, ub)
  value <- Upper(c(-0.1, 0, 0.1), H_k = H_k, dH_k = dH_k, X_k = X_k, Z_k = Z_k)
  expect_that( value[1] == value[2], is_true() )
  expect_that( value[2] == value[3], is_true() )
  print("Passed Upper function unit test")
})
#Test for lower
#The value of lower at 0 should equal to the value of h_function
g<- dnorm
test_that("Lower", {
  print("Testing Lower function now:")
  X_k <- c(-1,0,1)
  lb <- -1
  ub <- 1
  H_k <- h(X_k)
  dH_k <- grad(h,X_k)
  Z_k <- ConstructZ(X_k, h, lb, ub)
  value <- Lower(0,H_k = H_k, dH_k = dH_k, X_k = X_k)
  expect_that( value == h(0), is_true() )
  print("Passed Lower function unit test")
})
