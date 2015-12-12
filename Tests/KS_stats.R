
### This file generates Kolmogorov-Smirnoff statistics for evaluating the ARS algorithm for common log-concave distributions ###

### Exponential distribution ###
KS_exp = function(n,lambda){
  
  # generate samples from ARS function
  #x_fun = ars(dexp(x,lambda),n)
  
  # Use dummy exponential sampling function 
  x_fun = sample_exp(lambda,n)
  
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
KS_norm = function(n,mu,sigma){
  
  # generate samples from ARS function
  x_fun = ars(dnorm(x,mu,sigma),n)
  
  # construct empirical CDF
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
KS_unif = function(n,a,b){
  
  # generate samples from ARS function
  x_fun = ars(dunif(x,a,b),n)
  
  # construct empirical CDF
  x_fun = sort(x_fun)
  F_fun = (1:n)/n
  
  # Construct empirical CDF using built in R functions
  x_act = runif(n,mu,sigma)
  x_act = sort(x_act)
  
  F_act = sapply(x_act,function(x){ F_fun[which.min(abs(x-x_fun))] })
  
  
  # Compute Kolmogorov-Smirnoff statistic
  return(max(abs(F_fun - F_act)))
  
}


### Beta distribution ###
KS_beta = function(n,alpha,beta){
  
  # generate samples from ARS function
  x_fun = ars(dbeta(x,alpha,beta),n)
  
  # construct empirical CDF
  x_fun = sort(x_fun)
  F_fun = (1:n)/n
  
  # Construct empirical CDF using built in R functions
  x_act = runif(n,alpha,beta)
  x_act = sort(x_act)
  
  F_act = sapply(x_act,function(x){ F_fun[which.min(abs(x-x_fun))] })
  
  
  # Compute Kolmogorov-Smirnoff statistic
  return(max(abs(F_fun - F_act)))
  
}

### Gamma distribution ###
KS_gamma = function(n,shape,rate){
  
  # generate samples from ARS function
  x_fun = ars(dgamma(x,shape,rate),n)
  
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
KS_chi = function(n,df){
  
  # generate samples from ARS function
  x_fun = ars(dchisq(x,df),n)
  
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
KS_weibull = function(n,shape,scale){
  
  # generate samples from ARS function
  x_fun = ars(dweibull(x,shape,scale),n)
  
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
