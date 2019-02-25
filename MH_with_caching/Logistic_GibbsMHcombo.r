#Author: Arman Oganisian

library(LaplacesDemon)
library(invgamma)
library(MASS)
library(profvis)

################################################################################
### 0 - Simulate Data 
################################################################################
set.seed(10)
N<-10000

d<-data.frame(age_group=sample(x = c(0,1,2), size = N, replace = T))
d$age_1<-ifelse(d$age_group==1,1,0)
d$age_2<-ifelse(d$age_group==2,1,0)

d$trt<-rbinom(n = N, size = 1,prob = invlogit(0 + 2*d$age_1 + - 2*d$age_2))

d$y<-rbinom(n = N, size = 1,
            prob = invlogit(-1 + .7*d$age_1 + 1.1*d$age_2 + 1.1*d$trt))

X<-as.matrix(cbind(1,d[,2:4])) # model matrix
Y<-matrix(d$y, ncol=1) # outcome vector

p<-ncol(X)

################################################################################
### 1 - functions to sample from conditional posterior distributions
################################################################################

# unnormalized log posterior of beta vector
log_cond_post_beta<-function(beta, phi, lambda, X, Y){
  # calculate likelihood
  xb <- X %*% beta
  xb <- ifelse(xb>10, 10, ifelse(xb< (-10),-10, xb))
  p <- invlogit(xb)
  
  lik<- dbinom(x = Y, size = 1, prob = p, log = T)
  
  # calculate prior
  pr <- dnorm(x = beta, mean = lambda, sd = sqrt(phi), log = T)
  
  log_cond_post <- lik + pr
  return(log_cond_post)
}

# use Metropolis Hastings algorithm to sample from cond. post. of beta
mh_vanilla <- function(beta_0, phi, lambda, X, Y, mh_trials,jump_v){
  
  accept <- 0
  post_draws <- matrix(data = NA, nrow = mh_trials, ncol = length(beta_0))
  jump_cov <- jump_v*diag(p)
  
  for(i in 1:mh_trials){
    # draw from proposal distribution
    beta_c <- mvrnorm(1,beta_0,Sigma = jump_cov)
    
    # calculate ratio of conditional posterior densities
    r_num <- log_cond_post_beta(beta_c, phi, lambda, X, Y )
    r_denom <- log_cond_post_beta(beta_0, phi, lambda, X, Y )
    r <- exp(r_num - r_denom)
    rmin<-min(r,1)
    
    # accept or reject proposal
    if(rmin>=1){ 
      beta_0 <- beta_c 
      accept <- accept+1
    }else if(rbinom(1,1,rmin)==1){ 
      beta_0 <- beta_c
      accept <- accept+1
    }
    
    post_draws[i, ] <- beta_0
  }
  
  return(list(post_draws=post_draws, accept=accept/mh_trials)  )
}

mh_cache <- function(beta_0, phi, lambda, X, Y, mh_trials,jump_v){
  
  accept <- 0
  accept_flag <- 0
  post_draws <- matrix(data = NA, nrow = mh_trials, ncol = length(beta_0))
  
  eval_store <- matrix(data = NA, nrow = mh_trials, ncol = 2)
  jump_cov <- jump_v*diag(p)
  
  eval_curr <- log_cond_post_beta(beta_0, phi, lambda, X, Y )
  
  for(i in 1:mh_trials){
    # draw from proposal distribution
    beta_c <- mvrnorm(n = 1, mu = beta_0, Sigma = jump_cov)
    
    # calculate ratio of conditional posterior densities
    
    eval_prop <- log_cond_post_beta(beta_c, phi, lambda, X, Y )
    r <- exp(eval_prop - eval_curr)

    rmin<-min(r,1)
    
    # accept or reject proposal
    if(rmin>=1){ 
      beta_0 <- beta_c 
      accept <- accept+1
      eval_curr <- eval_prop
    }else if(rbinom(1,1,rmin)==1){ 
      beta_0 <- beta_c
      accept <- accept+1
      eval_curr <- eval_prop
    }
    
    post_draws[i, ] <- beta_0
  }
  
  return(list(post_draws=post_draws, accept=accept/mh_trials)  )
}


################################################################################
### 2 - Run Gibbs Sampler
################################################################################

# true hyperparameter values for betas
lambda<-c(0,0,0,0)
phi<-10


library(rbenchmark)

set.seed(1)
mh_draw <- mh_vanilla(beta_0 = c(0,0,0,0), 
                      phi = phi, lambda = lambda, 
                      X = X, Y = Y,
                      mh_trials=20000, jump_v=.2)

set.seed(1)
mh_draw <- mh_cache(beta_0 = c(0,0,0,0), 
                    phi = phi, lambda = lambda, 
                    X = X, Y = Y,
                    mh_trials=20000, jump_v=.2)


################################################################################
### 3 - Plot Results
################################################################################

par(mfrow=c(2,2))
plot(mh_draw$post_draws[,1],type='l',
     xlab='MCMC Iterations',
     ylab=c('Coefficient Draw'),
     main='Intercept')
abline(h=-1,col='red')
plot(mh_draw$post_draws[,2],type='l',
     xlab='MCMC Iterations',
     ylab=c('Coefficient Draw'),
     main='Age1')
abline(h=.7,col='red')
plot(mh_draw$post_draws[,3],type='l',
     xlab='MCMC Iterations',
     ylab=c('Coefficient Draw'),
     main='Age2')
abline(h=1.1,col='red')
plot(mh_draw$post_draws[,4],type='l',
     xlab='MCMC Iterations',
     ylab=c('Coefficient Draw'),
     main='Treatment')
abline(h=1.1,col='red')

