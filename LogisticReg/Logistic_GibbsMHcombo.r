#Author: Arman Oganisian

library(LaplacesDemon)
library(invgamma)
library(MASS)
library(profvis)

################################################################################
### 0 - Simulate Data 
################################################################################
set.seed(10)
N<-1000

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

### sample from conditional posterior of phi - conjugate!
rcond_post_phi<-function(beta, alpha, gamma, lambda, p){
  
  post_alpha<-alpha + p/2
  post_gamma<-gamma + .5*t(beta - lambda)%*%(beta - lambda)
  draw<-invgamma::rinvgamma(n = 1, shape = post_alpha, rate = post_gamma)
  
  return(draw)
}

# unnormalized log posterior of beta vector
log_cond_post_beta<-function(beta, phi, lambda, X, Y){
  # calculate likelihood
  lik<-0
  for(i in 1:length(Y)){
    xb <- X[i,] %*% beta
    xb<-ifelse(xb>10, 10, ifelse(xb< (-10),-10, xb))
    
    p_i<-invlogit(xb)
    
    lik<-lik + Y[i]*log(p_i) + (1 - Y[i])*log(1 - p_i)
  }

  # calculate prior 
  pr <- -.5 * (1/phi)*( t(beta - lambda)%*%(beta - lambda) )
  
  log_cond_post <- lik + pr
  return(log_cond_post)
}

# use Metropolis Hastings algorithm to sample from cond. post. of beta
rcond_post_beta_mh<-function(beta_0, phi, lambda, X, Y, mh_trials,jump_v){

  for(i in 1:mh_trials){
    # draw from proposal distribution
    beta_c <- mvrnorm(1,beta_0,Sigma = jump_v*diag(p))
    
    # calculate ratio of conditional posterior densities
    r_num <- log_cond_post_beta(beta_c, phi, lambda, X, Y )
    r_denom <- log_cond_post_beta(beta_0, phi, lambda, X, Y )
    r <- exp(r_num - r_denom)
    rmin<-min(r,1)
    
    # accept or reject proposal
    accept<-0
    if(rmin>=1){ 
      beta_0<-beta_c 
      accept<-1
    }else{ 
      if(rbinom(1,1,rmin)==1){ 
        beta_0<-beta_c
        accept<-1
      }
    }
  }
  
  return(c(new_beta=beta_0, accept=accept)  )
}

################################################################################
### 2 - Run Gibbs Sampler
################################################################################

### Gibbs Sampler
# true hyperparameter values for phi
alpha<-5
gamma<-2

# true hyperparameter values for betas
lambda<-c(0,0,0,0)
phi<-10000 # initialize 

# shell for storing results
gibbs_iter<-2000 + 1
gibbs_res<-matrix(nrow=gibbs_iter, ncol=p+2)

# initialize 
gibbs_res[1,1:p]<-c(0,0,0,0)

profvis(expr = {
for(i in 2:gibbs_iter){
  # sample from posterior of phi
  gibbs_res[i,p+1] <- rcond_post_phi(gibbs_res[i-1,1:p], 
                                     alpha, gamma, lambda, p)
  # sample from posterior of beta vector ( using MH )
  mh_draw <- rcond_post_beta_mh(gibbs_res[i-1,1:p], gibbs_res[i,p+1], 
                                lambda, X, Y, mh_trials=5, jump_v=.01)
  
  # store results
  gibbs_res[i,1:p] <- mh_draw[1:p]
  gibbs_res[i,p+2] <- mh_draw[p+1]
}
})

################################################################################
### 3 - Plot Results
################################################################################

par(mfrow=c(2,2))
plot(gibbs_res[,1],type='l',xlab='MCMC Iterations',ylab=c('Coefficient Draw'),
     main='Intercept')
abline(h=-1,col='red')
plot(gibbs_res[,2],type='l',xlab='MCMC Iterations',ylab=c('Coefficient Draw'),
     main='Age1')
abline(h=.7,col='red')
plot(gibbs_res[,3],type='l',xlab='MCMC Iterations',ylab=c('Coefficient Draw'),
     main='Age2')
abline(h=1.1,col='red')
plot(gibbs_res[,4],type='l',xlab='MCMC Iterations',ylab=c('Coefficient Draw'),
     main='Treatment')
abline(h=1.1,col='red')

# calculate posterior means and credible intervals
post_burn_trim<-gibbs_res[seq(1000,gibbs_iter,100),]
colMeans(post_burn_trim)
apply(post_burn_trim, 2, quantile, p=c(.025,.975))
