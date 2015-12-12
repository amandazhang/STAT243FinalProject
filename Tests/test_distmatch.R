#source("ars.R")


test_that("Empirical distributions match", {
  
  n = 10000
  
  print("Comparing exponential distributions:")
  expect_that( KS_exp(n,2), equals(0, tolerance  = 0.02))
  expect_that( KS_exp(n,5), equals(0, tolerance  = 0.02))
  expect_that( KS_exp(n,10), equals(0, tolerance  = 0.02))
  
  # placeholder to prevent testing, remove after ARS function completed
  if(1 == 0){
  
  print("Comparing normal distributions:")
  expect_that( KS_norm(n,0,1), equals(0, tolerance  = 0.02))
  expect_that( KS_norm(n,-5,1), equals(0, tolerance  = 0.02))
  expect_that( KS_norm(n,10,5), equals(0, tolerance  = 0.02))
  
  print("Comparing uniform distributions:")
  expect_that( KS_unif(n,0,1), equals(0, tolerance  = 0.02))
  expect_that( KS_unif(n,-10,10), equals(0, tolerance  = 0.02))
  expect_that( KS_unif(n,90,100), equals(0, tolerance  = 0.02))
  
  print("Comparing beta distributions:")
  expect_that( KS_beta(n,2,2), equals(0, tolerance  = 0.02))
  expect_that( KS_beta(n,2,5), equals(0, tolerance  = 0.02))
  expect_that( KS_beta(n,5,10), equals(0, tolerance  = 0.02))
  
  print("Comparing gamma distributions:")
  expect_that( KS_gamma(n,2,2), equals(0, tolerance  = 0.02))
  expect_that( KS_gamma(n,2,5), equals(0, tolerance  = 0.02))
  expect_that( KS_gamma(n,5,10), equals(0, tolerance  = 0.02))
  
  print("Comparing chi-squared distributions:")
  expect_that( KS_chi(n,3), equals(0, tolerance  = 0.02))
  expect_that( KS_chi(n,5), equals(0, tolerance  = 0.02))
  expect_that( KS_chi(n,10), equals(0, tolerance  = 0.02))
  
  print("Comparing weibull distributions:")
  expect_that( KS_weibull(n,1,2), equals(0, tolerance  = 0.02))
  expect_that( KS_weibull(n,2,5), equals(0, tolerance  = 0.02))
  expect_that( KS_weibull(n,10,2), equals(0, tolerance  = 0.02))
  
  }
  
})