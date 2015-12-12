#STAT 243 final project

# Initial value: number of starting points

# D: The domain of the function f(x)

#243 Project
library(numDeriv)

#g_x is the user defined function we are given
#T_k is the vector that contains k abscissae in D: x1,x2,...,xk
g <- dnorm
T_k <- c(0,1,2,3)
# lower bound
lb <- 0
# upper bound
ub <- 3

#log-function h(x)
h <- function(x){
  return(log(g(x)))
}
h_value <- h(T_k)

#Calculating z_j
z_function <- function(x = T_k, h){
  # x value is a vector which has length k-1. (1,2,3,...,k-1)
  # then in order to compute x_j - x_(j-1), create two vectors and subtract them.
  x_1 <- c(x, 0) #x1, x2, x3, ..., x_(k-1),0
  x_2 <- c(0, x) #0, x1, x2, ..., x_(k-2), x_(k-1)
  numerator <- h(x_1)-h(x_2)-x_1*grad(h,x_1)+x_2*grad(h,x_2)
  denominator <- grad(h,x_2)-grad(h,x_1)
  z <- numerator/denominator # z0, z1, ..., z_(k-1), z_k.
  # Set Z0 to be lower bound
  z[1] <- lb
  # Set Z_k to be upper bound
  z[length(x)+1] <- ub
  return(z)
}

z_value <- z_function(T_k, g)

# Upper hull formed from the tangents to h(x)
# u_k function
u_k_function <- function(x){
  #since the index of z starts at z0, up to zk, the length is actually k+1. So to match the index of z and x, we need to subtract 1 from index of z to get index of x respectively.
  j <- min(which(x <= z_value)-1) 
  return(h(T_k[j]) + (x - T_k[j]) * grad(h, T_k[j]))
}

# Lower hull:
# l_k function
l_k_function <- function(x){
  j <- min(which(x <= T_k))
  numerator <- (T_k[j]-x)*h(T_k[j-1]) + (x-T_k[j-1])*h(T_k[j])
  denominator <- T_k[j] - T_k[j-1]
  return(numerator/denominator)
}

#s_k part
exponential_u_k_function <- function(x){
  return(exp(u_k_function(x)))
}
c <- integrate(exponential_u_k_function, lower = lb, upper = ub)$value
s_k <- function(x){
  return(exponential_u_k_function(x)/c)
}
#inverse cdf of s_k
s_k_inverse_cdf <- function(x){
  j <- min(which(x < z_value))
  part_1 <- 1/grad(h,T_k[j])
  upper <- c*grad(h,T_k[j])
  lower <- exp(h(T_k[j])-T_k[j]*grad(h,T_k[j]))
  part_2 <- log((upper*x/lower)+exp(z_value[j-1]*grad(h,T_k[j])))
  return(part_1 * part_2)
}

upper <- c*grad(h,T_k[j])
a <- s_k_inverse_cdf(runif(1000)/c)