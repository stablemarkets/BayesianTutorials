
## Helper Functions 
### compute constant of tempered distribution
calc_k <- function(ftemp, temp){
  k <- integrate(ftemp, lower=-Inf, upper=Inf, temp=temp)
  return(k$value)
}

vanilla_mh <- function(log_target, theta_curr){
  theta_star <- rnorm(1, theta_curr, 1 )
  
  log_accept_prob <- log_target(theta_star) - log_target(theta_curr)
  accept <- log(runif(1)) < log_accept_prob
  if(is.na(accept)) browser()
  if(accept){ return(theta_star) }else{ return(theta_curr)}
}

### tempered distribution for various temperature, temp.
ftemp <- function(x, temp) ( (1/3)*dnorm(x,-20,1) + (1/3)*dnorm(x,0,1) + (1/3)*dnorm(x,20,1))^(1/temp)
k <- calc_k(ftemp, 100)

### log target distribution (ftemp with temp=1: ie untempered )
log_target <- function(x) log(ftemp(x, temp=1))
curve(ftemp(x,temp=1), from = -50, to=50, ylim=c(0,1), add=T)

k <- calc_k(ftemp, 100)
curve( (1/k)*ftemp(x,temp=100), from = -30, to=30, add=T, col='red')

k <- calc_k(ftemp, 100)
curve( (1/k)*ftemp(x,temp=100), from = -50, to=50, col='red')

#### Run Vanilla Metropolis Hastings

iter <- 10000
theta_shell <- numeric(iter)
theta_shell[1] <- 1

for(i in 2:iter){
  theta_shell[i] <- vanilla_mh(log_target = log_target, 
                               theta_curr = theta_shell[i-1])
}

par(mfrow=c(1,2))
plot(theta_shell, type='l')
hist(theta_shell, breaks=100, freq=F)
curve(ftemp(x,temp=1), from = -40, to=40,ylim=c(0,1), add=T, col='red', n = 100000,lwd=2)

## Parallel Tempering 

iter <- 10000
tempv <- c(1,200)

n_temps <- length(tempv)
temp_indx <- 1:n_temps

theta_shell <- matrix(0, nrow=iter, ncol=n_temps)
swap_shell <- matrix(nrow=iter, ncol=2)

for( i in 2:iter){
  
  ## update chains (potentially in parallel )
  for(t in temp_indx){
    log_target <- function(x) log( ftemp(x, temp=tempv[t]) )
    
    theta_shell[i, t] <-  vanilla_mh(log_target = log_target, 
                                     theta_curr = theta_shell[i-1, t])
  }
  
  ## propose swap, from swap_idx[1] (chain j) to swap_idx[2] (chain k)
  swap_idx <- sample(temp_indx, 2, replace = F)
  cj <- swap_idx[1]
  ck <- swap_idx[2]
  theta_j <- theta_shell[ i , cj]
  theta_k <- theta_shell[ i , ck]
  
  
  f1 <- tempv[cj]*( log_target( theta_j ) - log_target( theta_k ) )
  f2 <- tempv[ck]*( log_target( theta_k ) -log_target( theta_j )  )
  
  accept_prob <- min( c(1, exp(f1 + f2) ) )
    
  if( rbinom(1,1, accept_prob)==1 ){
    
    ## make the swap
    theta_shell[i, cj] <- theta_k
    theta_shell[i, ck] <- theta_j
    
    ## record the swap
    swap_shell[i, 1] <- cj
    swap_shell[i, 2] <- ck
    
  }
  
}

par(mfrow=c(1,2))
plot(theta_shell[,2], type='l', col='gray')
lines(theta_shell[,1], col='black')

hist(theta_shell[,1], breaks=100,freq = F)
curve(ftemp(x,temp=1), from = -40, to=40,ylim=c(0,1), add=T, col='red', n = 100000,lwd=2)


plot(theta_shell[,2], type='l')

sm <- swap_shell[which(swap_shell[,1]==1 |swap_shell[,2]==1) ,]

par(mfrow=c(1,1))
plot(theta_shell[,1], type='l')
points(1:iter, theta_shell[,1], pch=20, 
       col=ifelse( (!is.na(swap_shell[,1]) & !is.na(swap_shell[,2])  ) & (swap_shell[,1]==1 | swap_shell[,2]==1), 'red', NA ) )


hist(theta_shell[,2],breaks=100, freq=F)
k <- calc_k(ftemp, 100)
curve((1/k)*ftemp(x,temp=100), from = -60, to=40,ylim=c(0,1), add=T, col='red', n = 100000,lwd=2)

par(mfrow=c(1,2))

plot(theta_shell[s_shell==100], col='gray',
     ylab='MCMC Draws', xlab='Iteration', ylim=c(-40,60), pch=20,
     main='MCMC Chains')
points(theta_shell[s_shell==1], col='blue', pch=20)
legend('top', horiz = F, legend=c('Draws from Target Distribution', 
                                  'Draws from Tempered Distribution'), 
       col=c('blue', 'gray'), pch=c(20,20))


hist(theta_shell[s_shell==100], col='gray', breaks=100, freq=F, 
     ylim=c(0,.3), xlim=c(-20,20), xlab='MCMC Draws', main='MCMC Draws')
hist(theta_shell[s_shell==1], col='blue', breaks=100,add=T, freq=F)
legend('top', horiz = F, legend=c('Draws from Target Distribution', 
                                  'Draws from Tempered Distribution',
                                  'Target Distribution'),
       bty='n',
       col=c('blue', 'gray','red'), lty=c(NA,NA, 1),pch=c(15,15,NA))
curve(ftemp(x,temp=1), from = -20, to=20, ylim=c(0,1), add=T, col='red', n = 100000,lwd=2)

