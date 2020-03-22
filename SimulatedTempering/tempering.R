### ------------------------- Helper Functions  -----------------------------###

vanilla_mh <- function(log_target, theta_curr, prop_sd){
  theta_star <- rnorm(1, theta_curr, prop_sd )
  
  log_accept_prob <- log_target(theta_star) - log_target(theta_curr)
  accept <- log(runif(1)) < log_accept_prob
  if(is.na(accept)) browser()
  if(accept){ return(theta_star) }else{ return(theta_curr)}
}

### tempered distribution for various temperature, temp.
ftemp <- function(x, temp) ( (1/3)*dnorm(x,-20,1) + (1/3)*dnorm(x,0,1) + (1/3)*dnorm(x,20,1))^(1/temp)
log_target <- function(x) log(ftemp(x, temp=1))

### ------------------------- Parallel Tempering  ---------------------------###
set.seed(10)

iter <- 200
tempv <- c(1,100)

n_temps <- length(tempv)
temp_indx <- 1:n_temps

theta_shell <- matrix(0, nrow=iter, ncol=n_temps)
theta_shell_no_switch = matrix(0, nrow=iter, ncol=n_temps)
swap_shell <- numeric(length = iter)

for( i in 2:iter){
  
  ## update chains (potentially in parallel )
  for(t in temp_indx){
    log_target <- function(x) log( ftemp(x, temp=tempv[t]) )
    
    prop_sd = ifelse(t==1, 1, 10)
    
    theta_shell[i, t] <-  vanilla_mh(log_target = log_target, 
                                     theta_curr = theta_shell[i-1, t], 
                                     prop_sd = prop_sd)
    
  }
  theta_shell_no_switch[i, ] = theta_shell[i, ]
  
  ## propose swap, from swap_idx[1] (chain j) to swap_idx[2] (chain k)
  swap_idx <- sample(temp_indx, 2, replace = T)
  cj <- swap_idx[1]
  ck <- swap_idx[2]
  theta_j <- theta_shell[ i , cj]
  theta_k <- theta_shell[ i , ck]
  
  log_target <- function(x) log(ftemp(x, temp=1))
  f1 <- tempv[cj]*( log_target( theta_j ) - log_target( theta_k ) )
  f2 <- tempv[ck]*( log_target( theta_k ) - log_target( theta_j )  )
  
  accept_prob <- min( c(1, exp(f1 + f2) ) )
    
  if( rbinom(1,1, accept_prob)==1 ){
    
    ## make the swap
    theta_shell[i, cj] <- theta_k
    theta_shell[i, ck] <- theta_j
    
    ## record the swap
    swap_shell[i] = ifelse(cj!=ck, 1, 0)

  }
  
}

par(mfrow=c(1,2))
plot(theta_shell_no_switch[,2], type='l', col='gray', ylim=c(-40,40))
lines(theta_shell[,1], col='steelblue')
#lines(theta_shell_no_switch[,1], col='green')
points(1:iter, theta_shell[,1], pch=20, 
       col=ifelse(swap_shell==1, 'red', NA) )
lines(theta_shell_no_switch[,2], col='gray')

hist(theta_shell[,1], breaks=30,freq = F, 
     xlim=c(-25,25))
curve(ftemp(x,temp=1), from = -40, to=40,ylim=c(0,1), add=T, col='pink', n = 100000,lwd=2,
      ylim=c(0,.3))
