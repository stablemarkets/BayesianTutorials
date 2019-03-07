library(survival)
library(truncdist)
source("FunctionSourceCode.r")
################################################################################
### 0 - Simulate Data
################################################################################
set.seed(1)

n <- 1000
x1 <- rnorm(n)
x2 <- rnorm(n)
A <- rbinom(n, 1, .5)

X <- model.matrix(~ x1 + x2 + A)

true_beta <- matrix(c(1, .5, -.5, 1), ncol=1)
true_shape <- exp(X %*% true_beta)
true_scale <- 1

# simulate censoring and survival times
survt = rweibull(n, true_shape, true_scale) 
cent = rweibull(n, true_shape, true_scale)

## observed data:
#censoring indicator
delta <- cent < survt
survt[delta==1] <- cent[delta==1] # censor survival time.

## plot Kaplan-Meier Estimate
par(mfrow=c(1,1))
plot(survfit(Surv(survt, 1-delta) ~ 1))

################################################################################
### 1 - Run Augmented Sampler Accounting for Censoring 
################################################################################
iter <- 10000 # number of gibbs iterations

# shells for storing parameters
bshell <- matrix(NA, nrow=4, ncol = iter)
scaleshell <- numeric(iter)

scaleshell[1] <- c(0)
bshell[,1] <- c(0,0,0,0)

survt_all <- survt
n_miss <- sum(delta)
row_miss <- c(1:n)[delta]

par(mfrow=c(1,2))
plot(survfit(Surv(survt, 1-delta) ~ 1),conf.int = F, col='green',
     xlab=c('Time'),ylab='Survival Probability',
     main = 'Data augmentation with all subjects')

prop_covar <- diag(c(.0001,.0001,.0001,.0001))

for(i in 2:iter){
  ## sample from posterior of parameters, 
  ## conditional on observed and missing survival times
  bshell[,i] <- metrop_hastings(x_0 = bshell[,i-1, drop=F], iter = 1,
                                log_post_density = log_post_a,
                                prop_sigma = prop_covar, 
                                X=X, survt=survt_all, b=scaleshell[i-1] )$x_0
  
  scaleshell[i] <- metrop_hastings(x_0 = scaleshell[i-1], iter = 1,
                                   log_post_density = log_post_b,
                                   prop_sigma = matrix(.0001), 
                                   X=X, survt=survt_all, 
                                   beta=bshell[, i, drop=F] )$x_0
  
  ## sample from conditional posterior of missing survival times
  for(m in row_miss){
    survt_all[m] <- rtrunc(1, spec = 'weibull', 
                           a = survt[m], 
                           shape = exp( X[m,,drop=F] %*% bshell[,i, drop=F] ) , 
                           scale =  exp(scaleshell[i]) )
  }
  
  if(i>9500){
    post_draw <- rweibull(n, 
                          shape =  exp(X %*% bshell[,i,drop=F]), 
                          scale = rep(exp(scaleshell[i]), n)  )
    post_ecdf <- ecdf(post_draw)
    curve(1-post_ecdf(x), add=T, from=0, to=4, col='gray')
  }
  
}

prior_draw <- rweibull(1000, 
                       shape = exp(X %*% true_beta), 
                       scale = true_scale  )
prior_ecdf <- ecdf(prior_draw)
curve(1-prior_ecdf(x), add=T, from=0, to=4, col='red', lwd=2)

lines(survfit(Surv(survt, 1-delta) ~ 1),conf.int = F, col='green', lwd=2)
legend('topright',
       legend=c('Kaplan-Meier', 'Posterior Survival Draws', 'True Survival Curve'),
       col=c('green','gray','red'), bty='n', lty=c(1,1,1))

# par(mfrow=c(1,1))
# plot(exp(scaleshell[1000:iter]), type='l')
# abline(h=true_scale, col='red')
# 
# par(mfrow=c(2,2))
# plot(bshell[1,1000:iter], type='l')
# abline(h=true_beta[1,1], col='red')
# plot(bshell[2,1000:iter], type='l')
# abline(h=true_beta[2,1], col='red')
# plot(bshell[3,1000:iter], type='l')
# abline(h=true_beta[3,1], col='red')
# plot(bshell[4,1000:iter], type='l')
# abline(h=true_beta[4,1], col='red')


################################################################################
### 2 - Run Sampler with only uncensored patients
################################################################################
survt_obs <- survt[delta!=1]
X_obs <- X[delta!=1, ]

plot(survfit(Surv(survt, delta) ~ 1), conf.int=F,col='green',
     xlab=c('Time'),ylab='Survival Probability',
     main = 'Metropolis with only uncensored subjects')

for(i in 2:iter){
  bshell[,i] <- metrop_hastings(x_0 = bshell[,i-1, drop=F], iter = 1,
                                log_post_density = log_post_a,
                                prop_sigma = prop_covar, 
                                X=X_obs, survt=survt_obs, b=scaleshell[i-1] )$x_0
  
  scaleshell[i] <- metrop_hastings(x_0 = scaleshell[i-1], iter = 1,
                                    log_post_density = log_post_b,
                                    prop_sigma = matrix(.0001), 
                                    X=X_obs, survt=survt_obs, 
                                    beta=bshell[, i, drop=F] )$x_0
  
  if(i>9500){
    post_draw <- rweibull(nrow(X_obs), 
                          shape =  exp(X_obs %*% bshell[,i,drop=F]), 
                          scale = rep(exp(scaleshell[i]), n)  )
    post_ecdf <- ecdf(post_draw)
    curve(1-post_ecdf(x), add=T, from=0, to=4, col='gray')
  }
  
}

prior_draw <- rweibull(1000, 
                       shape = exp(X_obs %*% true_beta), 
                       scale = true_scale  )
prior_ecdf <- ecdf(prior_draw)
curve(1-prior_ecdf(x), add=T, from=0, to=4, col='red')

legend('topright',
       legend=c('Kaplan-Meier', 'Posterior Survival Draws', 'True Survival Curve'),
       col=c('green','gray','red'), bty='n', lty=c(1,1,1))

lines(survfit(Surv(survt, delta) ~ 1), conf.int=F,col='green')



