simulate_uni_hawkes <- function(mu, alpha, beta, n){
  arrivals <- numeric()
  eps <- 1e-6

  # recurrence function
  S <- function(k){
    if(k==1){
      return(1)
    }
    return( 1 + S(k-1)* exp(-beta * (arrivals[k]- arrivals[k-1])) )
  }

 # the function to be solved to obtain, u
 fu <- function(u, U, k){
   return(log(U) + mu*(u - arrivals[k]) + alpha/beta * S(k) * (1 - exp(-beta * (u - arrivals[k]))))
 }

 fu_prime <- function(u, U , k){
   return( mu + alpha * S(k) * (exp(-beta * (u - arrivals[k]))))
 }

 # iterative procedure to solve f(u)
 solve_u <- function(U,k){
   u_prev <- arrivals[k] - log(U)/mu
   u_next <- u_prev - fu(u_prev, U, k)/ fu_prime(u_prev, U, k)
   while(abs(u_next - u_prev) > eps){
     u_prev <- u_next
     u_next <- u_prev - fu(u_prev, U, k)/fu_prime(u_prev, U, k)
   }
   return(0.5 * (u_prev + u_next))
 }

 set.seed(1)
 t1 <- -log(runif(1))/mu
 arrivals <- c( arrivals, t1)
 k <- length(arrivals)
 while(k < n){
  U <- runif(1)
  t_next <- solve_u( U, k)
  arrivals <- c(arrivals, t_next)
  k <- length(arrivals)
 }
 return(arrivals)
}
