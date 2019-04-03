# sources:
# https://cran.r-project.org/web/packages/SurvRegCensCov/vignettes/weibull.pdf

library(survival)
library(truncdist)

setwd("/Users/aoganisi/Box Sync/Research/Analyses/BNP_CE/code/")
source("FunctionSourceCode.r")
################################################################################
### 0 - Simulate Data
################################################################################
set.seed(1)

n <- 1000
A <- rbinom(n, 1, .5)

X <- model.matrix(~ A)

true_beta <- (1/2)*matrix(c(-1/3, 2), ncol=1)
true_mu <- X %*% true_beta

true_sigma <- 1

true_alpha <- 1/true_sigma
true_lambda <- exp(-1*true_mu*true_alpha)

hist(rweibull(n, shape=true_alpha, scale = true_lambda), breaks=100)

# simulate censoring and survival times
survt = rweibull(n, shape=true_alpha, scale = true_lambda) 
cent = rweibull(n, shape=true_alpha, scale = true_lambda)

## observed data:
#censoring indicator
delta <- cent < survt
survt[delta==1] <- cent[delta==1] # censor survival time.

# survt_all will combine observed and imputed survival times.
survt_all <- survt

# count number of missing/censored survival times
n_miss <- sum(delta)
row_miss <- c(1:n)[delta] # index for which rows are censored

################################################################################
### 1 - Run Augmented Sampler Accounting for Censoring 
################################################################################
iter <- 10000 # number of gibbs iterations
burnin <- 9000 # burn-in iterations

# shells for storing parameters
hazard_ratio <- numeric(iter - burnin)

# initial values
beta_shell <- matrix(c(0,0), ncol=1)
lalpha_shell <- c(0)

prop_covar <- diag(c(.01,.01))

# plot stratified Kaplan-Meier
par(mfrow=c(1,1))
plot(survfit(Surv(survt, 1-delta) ~ A),conf.int = F, col=c('blue','red'),
             xlab=c('Time'),ylab='Survival Probability',
             main = 'Data augmentation with all subjects')

for(i in 2:iter){
  ## sample from posterior of parameters, 
  ## conditional on observed and missing survival times
  
  # metrop_hastings() is a custom function for generating a draw 
  # from conditional posterior of beta: log_post_beta
  beta_shell <- metrop_hastings(x_0 = beta_shell, 
                                    iter = 1,
                                    log_post_density = log_post_beta,
                                    prop_sigma = prop_covar, 
                                    X=X, survt=survt_all, 
                                    log_alpha=lalpha_shell )$x_0
  
  # sample from conditional posterior of alpha: log_post_alpha
  lalpha_shell <- metrop_hastings(x_0 = lalpha_shell, 
                                     iter = 1,
                                     log_post_density = log_post_alpha,
                                     prop_sigma = matrix(.001), 
                                     X=X, survt=survt_all, 
                                     beta=beta_shell)$x_0
  
  ## sample from conditional posterior of missing survival times
  mu_curr <-  X %*% beta_shell
  alpha_curr <- exp(lalpha_shell)
  
  for(m in row_miss){
    lambda_curr <- exp(-1*mu_curr[m]*alpha_curr)
    
    survt_all[m] <- rtrunc(1, spec = 'weibull', 
                           a = survt[m], 
                           shape =  alpha_curr, 
                           scale =  lambda_curr)
  }
  
  if(i>burnin){
    # plot 500 posterior survival curve draws for treated and placebo
    mu_trt <-  sum(beta_shell)
    mu_pbo <-  beta_shell[1]
    
    post_draw <- rweibull(n, shape = alpha_curr, scale = exp(-1*mu_trt*alpha_curr)  )
    post_ecdf <- ecdf(post_draw)
    curve(1-post_ecdf(x), add=T, from=0, to=15, col='lightblue')
    
    post_draw <- rweibull(n, shape = alpha_curr, scale = exp(-1*mu_pbo*alpha_curr)  )
    post_ecdf <- ecdf(post_draw)
    curve(1-post_ecdf(x), add=T, from=0, to=15, col='lightgray')
    
    # store hazard ratio
    hazard_ratio[i-burnin] <- exp(-beta_shell[2]*alpha_curr)
  }
  
}

# overlay KM curve and plot legend
lines(survfit(Surv(survt, 1-delta) ~ A),conf.int = T, col=c('black','blue'))
legend('topright', 
       legend = c('KM Curve and Intervals (TRT)',
                  'Posterior Survival Draws (TRT)',
                  'KM Curve and Intervals (PBO)',
                  'Posterior Survival Draws (PBO)'),
       col=c('black','gray','blue','lightblue'), 
       lty=c(1,0,1,0), pch=c(NA,15,NA,15), bty='n')

plot(hazard_ratio, type='l')
abline(h=exp(-true_beta[2]*true_alpha), col='red')

hist(hazard_ratio)
abline(v=exp(-true_beta[2]*true_alpha), col='red')

