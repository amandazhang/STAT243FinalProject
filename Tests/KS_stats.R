### This file generates Kolmogorov-Smirnoff statistics for evaluating the ARS algorithm for common log-concave distributions ###

### Exponential distribution ###
KS_exp = function(x_fun,lambda){

  # construct function
  g <- function(x) dexp(x,lambda)

  # generate samples from ARS function
  n = length(x_fun)

  # construct empirical CDF
  x_fun = sort(x_fun)
  F_fun = (1:n)/n

  # Construct empirical CDF using built in R functions
  x_act = rexp(n,lambda)
  x_act = sort(x_act)

  F_act = sapply(x_act,function(x){ F_fun[which.min(abs(x-x_fun))] })


  # Compute Kolmogorov-Smirnoff statistic
  return(max(abs(F_fun - F_act)))

}


### Normal distribution ###
KS_norm = function(x_fun,mu,sigma){

  # construct empirical CDF
  n = length(x_fun)
  x_fun = sort(x_fun)
  F_fun = (1:n)/n

  # Construct empirical CDF using built in R functions
  x_act = rnorm(n,mu,sigma)
  x_act = sort(x_act)

  F_act = sapply(x_act,function(x){ F_fun[which.min(abs(x-x_fun))] })


  # Compute Kolmogorov-Smirnoff statistic
  return(max(abs(F_fun - F_act)))

}


### Uniform distribution ###
KS_unif = function(x_fun,a,b){


  n = length(x_fun)

  # construct function
  g <- function(x) dunif(x,a,b)

  # generate samples from ARS function


  # construct empirical CDF
  x_fun = sort(x_fun)
  F_fun = (1:n)/n

  # Construct empirical CDF using built in R functions
  x_act = runif(n,a,b)
  x_act = sort(x_act)

  F_act = sapply(x_act,function(x){ F_fun[which.min(abs(x-x_fun))] })


  # Compute Kolmogorov-Smirnoff statistic
  return(max(abs(F_fun - F_act)))

}


### Beta distribution ###
KS_beta = function(x_fun,alpha,beta){

  # construct function
  g <- function(x) dbeta(x,alpha,beta)

  n <- length(x_fun)


  # construct empirical CDF
  x_fun = sort(x_fun)
  F_fun = (1:n)/n

  # Construct empirical CDF using built in R functions
  x_act = rbeta(n,alpha,beta)
  x_act = sort(x_act)

  F_act = sapply(x_act,function(x){ F_fun[which.min(abs(x-x_fun))] })


  # Compute Kolmogorov-Smirnoff statistic
  return(max(abs(F_fun - F_act)))

}

### Gamma distribution ###
KS_gamma = function(x_fun,shape,rate){


  # construct function
  g <- function(x) dgamma(x,shape,rate)

  n = length(x_fun)

  # construct empirical CDF
  x_fun = sort(x_fun)
  F_fun = (1:n)/n

  # Construct empirical CDF using built in R functions
  x_act = rgamma(n,shape,rate)
  x_act = sort(x_act)

  F_act = sapply(x_act,function(x){ F_fun[which.min(abs(x-x_fun))] })


  # Compute Kolmogorov-Smirnoff statistic
  return(max(abs(F_fun - F_act)))

}

### Chi-square distribution ###
KS_chi = function(x_fun,df){


  # construct function
  g <- function(x) dchisq(x,df)

  n = length(x_fun)

  # construct empirical CDF
  x_fun = sort(x_fun)
  F_fun = (1:n)/n

  # Construct empirical CDF using built in R functions
  x_act = rchisq(n,df)
  x_act = sort(x_act)

  F_act = sapply(x_act,function(x){ F_fun[which.min(abs(x-x_fun))] })


  # Compute Kolmogorov-Smirnoff statistic
  return(max(abs(F_fun - F_act)))

}

### Weibull distribution ###
KS_weibull = function(x_fun,shape,scale){

  # construct function
  g <- function(x) dweibull(x,shape,scale)

  # generate samples from ARS function
  n = length(x_fun)

  # construct empirical CDF
  x_fun = sort(x_fun)
  F_fun = (1:n)/n

  # Construct empirical CDF using built in R functions
  x_act = rweibull(n,shape,scale)
  x_act = sort(x_act)

  F_act = sapply(x_act,function(x){ F_fun[which.min(abs(x-x_fun))] })


  # Compute Kolmogorov-Smirnoff statistic
  return(max(abs(F_fun - F_act)))

}
