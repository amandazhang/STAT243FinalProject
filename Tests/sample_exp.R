sample_exp  <- function(lambda,n){
  
  
  y = runif(n)
  x = -log(1-y)/lambda
  
  return(x)
  
  
  
  
  
}