rcond_post_beta <- function(nparams, y, xm, psi, 
                            beta_prior_mean, beta_prior_var){
  
  mu_beta <- beta_prior_mean
  v_beta <- diag(beta_prior_var)
  v_beta_inv <- diag((1/beta_prior_var))
  
  xtx <- t(xm)%*%xm
    
  post_cov <- solve( v_beta_inv + (1/psi)*xtx)
  post_mean <- post_cov %*% (v_beta_inv %*% mu_beta + (1/psi)*t(xm)%*%y )
    
  draw <- rmvnorm(n = 1, mean = post_mean, sigma = post_cov)
  return(draw)  
}

sim_dat <- function(n=10000){
  x1 <- rbinom(n = n, size = 1, prob = .5)
  x2 <- rnorm(n = n, mean = 0, sd = 10)
  
  p <- invlogit(1 + -2*x1 + 1*x2)
  
  y <- rbinom(n = n, size = 1, prob = p)
  
  d <- data.frame(y=y, x1=x1, x2=x2)
  return(d)
}
