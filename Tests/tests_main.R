rm(list=ls())
source('KS_stats.R')
source('ARS.R')

test_that("Empirical distribution matches", {

  print("Testing normal distribution now:")
  g <- function(x) dnorm(x,0,1)
  x <- ars(g,5000,-10,10,batchsize=100)
  print(expect_that( KS_norm(x,0,1), equals(0, tolerance  = 0.05)))
  cat("\n\n")

  print("Testing gamma distribution now:")
  g <- function(x) dgamma(x,2,2)
  x <- ars(g,5000,0,10,batchsize=100)
  print(expect_that( KS_gamma(x,2,2), equals(0, tolerance  = 0.05)))
  cat("\n\n")

  print("Testing uniform distribution now:")
  g <- function(x) dunif(x,0,1)
  x <- ars(g,5000,0,1,batchsize=100)
  print(expect_that( KS_unif(x,0,1), equals(0, tolerance  = 0.05)))
  cat("\n\n")

  print("Testing beta distribution now:")
  g <- function(x) dbeta(x,2,2)
  x <- ars(g,5000,0.01,0.99,batchsize=100)
  print(expect_that( KS_beta(x,2,2), equals(0, tolerance  = 0.05)))
  cat("\n\n")

  print("Testing chi-square distribution now:")
  g <- function(x) dchisq(x,3)
  x <- ars(g,5000,0.01,10,batchsize=1000)
  print(expect_that( KS_chi(x,3), equals(0, tolerance  = 0.05)))


})
