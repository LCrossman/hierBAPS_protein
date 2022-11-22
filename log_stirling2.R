log_stirling2 <- function(n, k){
  if(!is.numeric(n)) stop("n is not numeric!")
  if(!is.numeric(k)) stop("k is not numeric!")
  if(k>n) stop("k must be less than n!")
  
  v <- n/k
  G <- lambertW(-v*exp(-v))
  
  lS2 <- log(sqrt((v-1)/(v*(1-G)))) +
    (n-k)*(log(v-1)-log(v-G)) +
    n*log(k)-k*log(n) +
    k*(1-G) +
    lchoose(n, k)
  return(lS2)
}
