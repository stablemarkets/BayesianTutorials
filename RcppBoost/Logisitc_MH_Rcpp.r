#Author: Arman Oganisian
library(LaplacesDemon)
library(invgamma)
library(MASS)
library(tidyr)
library(dplyr)
library(ggplot2)
library(microbenchmark)
library(Rcpp)

sourceCpp('log_post.cpp')

################################################################################
### 0 - Simulate Data 
################################################################################
set.seed(10)

sim_dat <- function(N){
  d<-data.frame(age_group=sample(x = c(0,1,2), size = N, replace = T))
  d$age_1<-ifelse(d$age_group==1,1,0)
  d$age_2<-ifelse(d$age_group==2,1,0)
  
  d$trt<-rbinom(n = N, size = 1,prob = invlogit(0 + 2*d$age_1 + - 2*d$age_2))
  
  d$y<-rbinom(n = N, size = 1,
              prob = invlogit(-1 + .7*d$age_1 + 1.1*d$age_2 + 1.1*d$trt))
  
  X<-as.matrix(cbind(1,d[,2:4])) # model matrix
  Y<-matrix(d$y, ncol=1) # outcome vector
  return(list(X=X, Y=Y))  
}

d <- sim_dat(N=1000)
X <- d$X
Y <- d$Y

################################################################################
### 1 - functions to sample from conditional posterior distributions
################################################################################

# unnormalized log posterior of beta vector
log_posterior<-function(beta, X, Y){

    # calculate likelihood
  xb <- X %*% beta
  xb <- ifelse(xb>10, 10, ifelse( xb< (-10) ,-10, xb))
  p_i <- invlogit(xb)
  
  lik <- sum(dbern(Y, p_i, log = T))
  
  # calculate prior 
  pr <- dmvn(x = beta, mu = rep(0,p), Sigma = (1000^2)*diag(p), log = T)
  
  log_cond_post <- lik + pr
  return(log_cond_post)
}

# Metropolis-Hastings Sampler using log_posterior(),
# which is the log posterior coded in R.
sample_mh<-function(X, Y, iter, jump_v){
  
  # create shells
  p <- ncol(X)
  beta_shell <- matrix(NA, nrow = iter, ncol = p)
  accept_shell <- numeric(length = iter)
  
  # starting values
  beta_shell[1,] <- rep(10, p)
  
  for(i in 2:iter){
    beta_0 <- beta_shell[i-1, ]
    
    # draw from proposal distribution
    beta_c <- mvrnorm(n = 1, beta_0, Sigma = jump_v*diag(p))
    
    # calculate ratio of conditional posterior densities
    r_num <- log_posterior(beta_c, X, Y )
    r_denom <- log_posterior(beta_0, X, Y )
    
    # calculate acceptance probability
    r <- exp(r_num - r_denom)
    rmin<-min(r,1)
    
    # accept or reject proposal
    if( rbinom(1,1,rmin) == 1 ){ 
      beta_shell[i, ] <- beta_c
    }else{
      beta_shell[i, ] <- beta_0
    }
    accept_shell[i] <- rmin
    
  }
  colnames(beta_shell) <- colnames(X)
  colnames(beta_shell)[1] <- 'intercept'
  return(list(beta_shell, accept_shell) )
}

# Metropolis-Hastings Sampler using log_post(),
# which is the log posterior coded in C++.

sample_mh_cpp <-function(X, Y, iter, jump_v){
  # create shells
  p <- ncol(X)
  beta_shell <- matrix(NA, nrow = iter, ncol = p)
  accept_shell <- numeric(length = iter)
  
  # starting values
  beta_shell[1,] <- rep(10, p)
  
  for(i in 2:iter){
    beta_0 <- beta_shell[i-1, ]
    
    # draw from proposal distribution
    beta_c <- mvrnorm(n = 1, beta_0, Sigma = jump_v*diag(p))
    
    # calculate ratio of conditional posterior densities
    r_num <- log_post(beta_c, Y, X )
    r_denom <- log_post(beta_0, Y, X )
    
    # calculate acceptance probability
    r <- exp(r_num - r_denom)
    rmin<-min(r,1)
    
    # accept or reject proposal
    if( rbinom(1,1,rmin) == 1 ){ 
      beta_shell[i, ] <- beta_c
    }else{
      beta_shell[i, ] <- beta_0
    }
    accept_shell[i] <- rmin
    
  }
  colnames(beta_shell) <- colnames(X)
  colnames(beta_shell)[1] <- 'intercept'
  return(list(beta_shell, accept_shell) )
}

################################################################################
### 2 - Test the Samplers
################################################################################

burnin <- 1000
iter <- 100000
p <- ncol(X)

res_mh_cpp <- sample_mh_cpp(X, Y, iter = iter, jump_v = .03)

par(mfrow=c(2,2))
plot(res_mh_cpp[[1]][burnin:iter,'intercept'], type='l',
     xlab='MH Iteration', ylab='Posterior Draw', main='Intercept')
abline(h= -1, col='red')
plot(res_mh_cpp[[1]][burnin:iter,'age_1'], type='l',
     xlab='MH Iteration', ylab='Posterior Draw', main='age1')
abline(h= .7, col='red')
plot(res_mh_cpp[[1]][burnin:iter,'age_2'], type='l',
     xlab='MH Iteration', ylab='Posterior Draw', main='age2')
abline(h= 1.1, col='red')
plot(res_mh_cpp[[1]][burnin:iter,'trt'], type='l',
     xlab='MH Iteration', ylab='Posterior Draw', main='trt')
abline(h= 1.1, col='red')

par(mfrow=c(1,1))
plot(cumsum(res_mh_cpp[[2]])/1:iter, type='l',
     xlab='MH Iteration', ylab='Cumulative Average Acceptance Rate', 
     main='Acceptance Rate Over Sampling Run')
abline(h= 1.1, col='red')

################################################################################
### 3 - Benchmarks
################################################################################
iter <- 10000

ss <- c(100, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000)

rel_time <-numeric(length = length(ss))
for(i in 1:length(ss) ){
  d <- sim_dat(N = ss[i])
  X <- d$X
  Y <- d$Y
  
  bench<-microbenchmark(R_MH = sample_mh(X, Y, iter = iter, jump_v = .03),
                        Cpp_MH = sample_mh_cpp(X, Y, iter = iter, jump_v = .03),
                        times = 2)
  bench_sum <- summary(bench)
  r_time <- bench_sum$mean[bench_sum$expr=='R_MH']
  rcpp_time <- bench_sum$mean[bench_sum$expr=='Cpp_MH']
  rel_time[i] <- r_time/rcpp_time
}

################################################################################
### 4 - Plot Results
################################################################################

plot(ss, rel_time, type='l')


