#243 Project
library(numDeriv)
#g is the user defined function we are given
g <- dnorm
x_value <- c(-1,0,1)
lb <- -2
ub <- 2

#log-function
h <- function(x){
  return(log(g(x)))
}
h_value <- h(x_value)
dh_value <- grad(h, x_value)
#z-values
z_function <- function(x, h){
  z <- c(mode = "numeric",length = length(x))
  x_1 <- c(x, 0)
  x_2 <- c(0, x)
  upper <- h(x_1)-h(x_2)-x_1*grad(h,x_1)+x_2*grad(h,x_2)
  lower <- grad(h,x_2)-grad(h,x_1)
  tmp <- upper/lower
  tmp[1] <- lb
  tmp[length(x)+1] <- ub
  return(tmp)
}

z_value <- z_function(x_value, h)

#u_k function, some modification required
u_k_function <- function(x){
  j <- min(which(x < z_value))
  return(h_value[j] + (x - x_value[j]) * dh_value[j])
}

#l_k part, some modification required
l_k_function <- function(x){
  j <- min(which(x < z_value))
  upper <- (x_value[j+1]-x)*h_value[j] + (x-x_value[j])*h_value[j+1]
  lower <- x_value[j+1] - x_value[j]
  return(upper/lower)
}

#s_k part
exponential_u_k_function <- function(x){
  return(exp(u_k_function(x)))
}

#calculate constant c
c_function <- function(x_value, z_value, h_value, dh_value){
  upper <- exp(h_value - x_value * dh_value)
  lower <- dh_value
  part_1 <- upper/lower
  part_2 <- exp(z_value[2:length(z_value)]*dh_value) - exp(z_value[1:length(z_value)-1]*dh_value)
  result <- part_1 * part_2
  j <- which(TRUE == is.nan(result))
  if(j != 0){
    result[j] <- (z_value[j] - z_value[j-1])*exp(h_value[j])
    return(result)
  }else{
    return(result)
  }
}

c <- c_function(x_value, z_value, h_value, dh_value)
constant <- sum(c)
cum <- cumsum(c)
cum <- c(0,cum)

s_k <- function(x){
  return(exponential_u_k_function(x)/constant)
}


#inverse CDF part
random_number_function <- function(ns){
  runif(ns)*cum[length(cum)]
}

random_number <- random_number_function(50)

sample_function <- function(u,cum,constant){
  j <- max(which(u > cum))
  if(dh_value[j] == 0){
    x_star <- runif(1, z_value[j],z_value[j+1])
    return(x_star)
  }else{
    uc <- u*constant
    upper_1 <-(uc - cum[j])*dh_value[j]
    lower_1 <- exp(h_value[j]) - x_value[j]*dh_value[j]
    part_1 <- upper_1/lower_1
    upper_2 <- exp(dh_value[j] * z_value[j])+part_1
    log_upper_2 <- log(upper_2)
    lower_2 <- dh_value[j]
    x_star <- log_upper_2/lower_2
    return(x_star)
  }
}
x_star <- sapply(random_number, sample_function, cum = cum, constant = constant)



















#Main function

#h_function
h <- function(x){
  return(log(g(x)))
}

#z_function
z_function <- function(x, h, lb, ub){
  z <- c(mode = "numeric",length = length(x))
  x_1 <- c(x, 0)
  x_2 <- c(0, x)
  upper <- h(x_1)-h(x_2)-x_1*grad(h,x_1)+x_2*grad(h,x_2)
  lower <- grad(h,x_2)-grad(h,x_1)
  tmp <- upper/lower
  tmp[1] <- lb
  tmp[length(x)+1] <- ub
  return(tmp)
}

#u_k function, some modification required
u_k_function <- function(x,x_value, dh_value, h_value, z_value){
  j <-max(which(x > z_value))
  return(h_value[j] + (x - x_value[j]) * dh_value[j])
}

#l_k part, some modification required
l_k_function <- function(x, x_value, dh_value, h_value, z_value){
  j <-max(which(x > z_value))
  if(j == 1 | j == length(x_value)){
    return(-Inf)
  }else{
    upper <- (x_value[j]-x)*h_value[j-1] + (x-x_value[j-1])*h_value[j]
    lower <- x_value[j] - x_value[j-1]
    return(upper/lower)
  }
}

#calculate constant c
c_function <- function(x_value, z_value, h_value, dh_value){
  upper <- exp(h_value - x_value * dh_value)
  lower <- dh_value
  part_1 <- upper/lower
  part_2 <- exp(z_value[2:length(z_value)]*dh_value) - exp(z_value[1:length(z_value)-1]*dh_value)
  result <- part_1 * part_2
  if(sum(is.nan(result))){
    j <- which(TRUE == is.nan(result))
    result[j] <- (z_value[j] - z_value[j-1])*exp(h_value[j])
    return(result)
  }else{
    return(result)
  }
}

#random_number-generator
random_number_function <- function(ns, cum){
  runif(ns)*cum[length(cum)]
}

#This is the sample function, which gives us the x_star
sample_function <- function(u,cum,constant, dh_value, z_value, x_value, h_value){
  j <- max(which(u > cum))
  if(dh_value[j] == 0){
    x_star <- runif(1, z_value[j],z_value[j+1])
    return(x_star)
  }else{
    uc <- u*constant
    upper_1 <-(uc - cum[j])*dh_value[j]
    lower_1 <- exp(h_value[j]) - x_value[j]*dh_value[j]
    part_1 <- upper_1/lower_1
    upper_2 <- exp(dh_value[j] * z_value[j])+part_1
    log_upper_2 <- log(upper_2)
    lower_2 <- dh_value[j]
    x_star <- log_upper_2/lower_2
    return(x_star)
  }
}

ars_function <- function(g, ns, x_values, lb = -2, ub = 2){
  #h function
  h <- function(x){
    return(log(g(x)))
  }
  sample <- vector(mode = "numeric", length = 0)
  #x_value
  x_value <- x_values
  while(length(sample) < ns){
    
    print(length(sample))
    #h_value
    h_value <- h(x_value)
    dh_value <- grad(h, x_value)
    #z_value
    z_value <- z_function(x_value, h, lb, ub)
    
    #c as a vector storing all different areas under u(x)
    c <- c_function(x_value, z_value, h_value, dh_value)
    #normalizing constant
    constant <- sum(c)
    #cummulative sum
    cum <- cumsum(c)
    #add a 0 into cum, for further index use
    cum <- c(0,cum)
    #This gets us some random numbers
    random_number <- random_number_function(ns,cum)
    #This gives us the x-stars
    x_star <- sapply(random_number, sample_function, cum = cum, constant = constant, dh_value = dh_value, h_value=h_value, z_value=z_value, x_value = x_value)
    
    #Here we generate some w from uniform(0,1)
    w <- runif(ns)
    #squeezing test
    l_k_value <- sapply(x_star, l_k_function, dh_value = dh_value, h_value=h_value, z_value=z_value, x_value = x_value)
    u_k_value <- sapply(x_star, u_k_function, dh_value = dh_value, h_value=h_value, z_value=z_value, x_value = x_value)
    w_right <- exp(l_k_value - u_k_value)
    position_reject <- which(w > w_right)
    if (min(position_reject) != 1){
      sample <- x_star[1:min(position_reject)-1]
    }
    check_point <- min(position_reject)
    #rejection test
    h_value_check <- h(x_star[check_point])
    u_k_value_2 <- u_k_function(x = x_star[check_point],x_value = x_value, dh_value = dh_value, h_value = h_value, z_value = z_value)
    if(w[check_point] <= exp(h_value_check)-u_k_value_2){
      sample <- c(sample, x_star[check_point])
    }else{
      x_value <- c(x_value, x_star[check_point])
    }
    
  }
}

a <- c(-1,0,1)
sample <- ars_function(dnorm, 100, a, -2 , 2)