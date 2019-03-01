## Author: Arman Oganisian

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
mh_vanilla <- function(beta_0, phi, lambda, X, Y, mh_trials,jump_v, p){
  
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

mh_cache <- function(beta_0, phi, lambda, X, Y, mh_trials,jump_v, p){
  
  accept <- 0
  post_draws <- matrix(data = NA, nrow = mh_trials, ncol = length(beta_0))
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