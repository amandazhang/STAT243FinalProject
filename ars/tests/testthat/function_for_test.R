#### Adaptive Rejection Sampling ####
library(numDeriv)
## Take log of input function ##
h <- function(x){
  return(log(g(x)))
}

## Function that tests for log concavity ##
logconc <- function(h,lb,ub) {

  #log concave test here
  p <- runif(10000, min = 0, max = 1)

  #tests log concavity
  result <- h(p*lb + (1-p)*ub) >= p*h(lb) + (1-p)*h(ub)

  #result is 10000 TRUE's for the 10000 in p
  t <- rep(TRUE, 10000)
  finalreturn <- all.equal(t, result)

  #if returns TRUE then all elements are equal and 10000 reps were log concave so log concave
  if (finalreturn == FALSE) {
    warning("Density input to ars likely not log-concave. Please check and try again.")
  }
}

## Test for uniform case ##
utest <- function(g, lb, ub) {
  d <- (ub-lb)/100
  delta <- (ub - lb)/100
  x <- seq(lb + delta, ub - delta, by = d)
  dgvals <- grad(g, x)
  t <- rep(0,length(x))
  finalreturn <- all.equal(t, dgvals)
  if (finalreturn == TRUE) {
    uniformcase <- TRUE
  } else {
    uniformcase <- FALSE
  }
  return(uniformcase)
}

## Find optimum function and X_init ##
findmax <- function(h,lb,ub) {
  delta <- (ub - lb)/1000
  max <- optimize(f = h, interval = c(lb, ub), lower = lb, upper = ub, maximum = TRUE)$maximum
  #taking care of exp case
  if (abs(max - ub) < .0001) {
    rp <- max
    mid <- max - .5*delta
    lp <- (max - delta)
    X_init <- c(lp,mid,rp)
  } else if (abs(max - lb) < .0001) {
    rp <- (max + delta)
    mid <- max + .5*delta
    lp <- max
    X_init <- c(lp,rp)
  } else {
    rp <- (max + delta)
    lp <- (max - delta)
    X_init <- c(lp,max,rp)
  }
  return(X_init)
}

## Compute intersection points of tangents, z ##
ConstructZ <- function(x,h,lb,ub){
  k = length(x)
  x_j1 <- c(0,x)
  x_j <- c(x,0)
  Num <- h(x_j1)-h(x_j)-x_j1*grad(h,x_j1)+x_j*grad(h,x_j)
  Den <- grad(h,x_j)-grad(h,x_j1)

  if(sum(abs(Den),na.rm = TRUE)<1e-6){

    tmp <-  (x[-1] + x[-k])/2
    tmp <- c(lb+1e-4,tmp,ub)

  }else{

    tmp <- Num/Den
    tmp[1] <- lb
    tmp[length(x)+1] <- ub
  }
  return(tmp)
}

## Lower bound, l(x) ##
Lower <- function(x,H_k,dH_k,X_k){

  if(x < min(X_k) | x > max(X_k)){
    return(-Inf)
  }else{
    j <- max(which(x >= X_k))
    Num <- (X_k[j+1]-x)*H_k[j] + (x-X_k[j])*H_k[j+1]
    Den <- X_k[j+1] - X_k[j]
    return(Num/Den)
  }
}

## Upper bound, u(x) ##
Upper <- function(x,H_k,dH_k,X_k,Z_k){
  j <- min(which(x < Z_k)-1)
  return(H_k[j] + (x - X_k[j]) * dH_k[j])
}

## Exponentiated upper bound, exp(u(x)) ##
ExpUpper <- function(x,H_k,dH_k,X_k,Z_k){
  return(exp(Upper(x,H_k,dH_k,X_k,Z_k)))
}

## Normalized version of exponentiated upper bound, s(x) ##
Envelope <- function(x,C,H_k,dH_k,X_k,Z_k){

  return(ExpUpper(x,H_k,dH_k,X_k,Z_k)/C)
}

## Generate candidate samples by sampling from envelope s(x) ##
# We sample from the inverse CDF of the envelope by first randomly selection a segment of the envelope, and then treating
# the selected segment as a PDF supported by the adjacent Z_k points
GenCandidates <- function(u,cum_area_env,H_k,X_k,dH_k,Z_k,areas_u){
  j <- max(which(u > cum_area_env))
  if(dH_k[j] == 0){
    x <- runif(1, Z_k[j],Z_k[j+1])
    return(x)
  }else{

    # Sample from uniform random
    w = runif(1)

    # Scale seed value w to area of the selected segment, since area under segment is not equal to 1
    w_sc = w*(1/areas_u[j])*exp(H_k[j] - X_k[j]*dH_k[j])*(exp(dH_k[j]*Z_k[j+1]) -exp(dH_k[j]*Z_k[j]))

    # Use inverse CDF of selected segment to generate a sample
    x = (1/dH_k[j])*log(w_sc*areas_u[j]/(exp(H_k[j] - X_k[j]*dH_k[j])) + exp(Z_k[j]*dH_k[j]))
  }
  return(x)
}

## Apply squeeze and rejection test to each candidate sample ##
RejectionTest <- function(x,H_k,dH_k,X_k,Z_k){

  # Generate random seed
  w = runif(1)

  # Initialize squeeze and reject tests, and indicator for adding point
  squeeze = FALSE
  accept = FALSE
  add = FALSE

  # Compute threshold values for failing squeeze and rejection tests
  l_threshold = exp(Lower(x,H_k,dH_k,X_k) - Upper(x,H_k,dH_k,X_k,Z_k))
  u_threshold = exp(h(x) - Upper(x,H_k,dH_k,X_k,Z_k))

  #print(l_threshold)
  if( w <= l_threshold){

    squeeze = TRUE
    accept = TRUE

  }else if(w <= u_threshold){

    squeeze = FALSE
    accept = TRUE

  }else{

    accept = FALSE
  }

  # Determine whether to add point to abscissae
  if(squeeze*accept==FALSE) add = TRUE

  # Return boolean indicating whether to accept candidate sample point
  return(list(rej=squeeze+accept,add=add))
}

#' Adaptive Rejection Sampling (ars) Function
#'
#' Adaptive rejection sampling function that creates a sample
#' from a density function based on the Gilks(1992) paper.
#'
#' @param g probability density function, as a function of x (e.g. g <- function(x) dnorm(x,1,2))
#' @param n number of values the final sample should contain
#' @param lb lower bound to evaluate the density on
#' @param ub upper bound to evaluate the density on
#' @param X_init optional vector of initial values, default = NULL
#' @param batchsize optional argument to adjust number of x* values to evaluate at once
#' @return sample of with specified (n) elements
#' @export

##########################################################################

### Main ARS function ###
ars <- function(g,n,lb,ub,X_init = NULL, batchsize = round(n/100)){

  # Test for uniform case
  uniftest <- utest(g, lb, ub)

  if(uniftest == TRUE) {
    x_all <- runif(n,lb,ub)
    return(x_all)
  }
  # Compute log of input function
  h <- function(x){
    return(log(g(x)))
  }

  if (is.null(X_init) == TRUE) {
    X_k <- findmax(h,lb,ub)
  } else {
    X_k <- X_init
  }

  # Initialize abscissae and sample points
  x_all = NULL

  while(length(x_all)<n){

    # Compute intersection points
    Z_k <- ConstructZ(X_k,h,lb,ub)

    # Store h and h' values
    H_k <- h(X_k)
    dH_k <- grad(h, X_k)

    # Calculate areas under exponential upper bound function for normalization purposes
    areas_u = unlist(sapply(2:length(Z_k),function(i){integrate(ExpUpper,Z_k[i-1],Z_k[i],H_k,dH_k,X_k,Z_k)})[1,])
    C = sum(areas_u)

    # compute cumulative areas under envelope, which will sum to 1
    areas_env = areas_u/C
    cum_area_env <- cumsum(areas_env)
    cum_area_env <- c(0,cum_area_env)

    # Generate seeds for Inverse CDF method
    seeds <-runif(batchsize)

    # Generate candidate samples
    x_candidates <- sapply(seeds, GenCandidates, cum_area_env = cum_area_env, H_k = H_k, X_k = X_k,dH_k = dH_k,Z_k = Z_k,areas_u = areas_u)

    # Apply squeeze and rejection tests
    test_flag <- sapply(x_candidates,RejectionTest, H_k = H_k, X_k = X_k, dH_k = dH_k,Z_k = Z_k)
    keep_sample_flag <- test_flag[1,]
    add_X_k_flag <- test_flag[2,]

    # Filter out rejected points and update full set of samples
    x_keep <- x_candidates[keep_sample_flag>0]
    x_all = c(x_keep,x_all)

    # Identify samples to add to X_k
    X_new_k = x_candidates[add_X_k_flag>0]
    X_k = sort(c(X_k,X_new_k))
  }
  return(x_all)
}
g <- dnorm
