library(rstan)
setwd("/Users/aoganisi/Dropbox/Stable Markets/BayesianTutorials/PartialPooling")

################################################################################
#####                       Simulate Data                                  #####
################################################################################

d_A<-rbinom(n = 5, size = 1, prob =  .2)
d_B<-rbinom(n = 10, size = 1, prob = .3)
d_C<-rbinom(n = 20, size = 1, prob = .5)
d_D<-rbinom(n = 30, size = 1, prob = .7)
d_E<-rbinom(n = 40, size = 1, prob = .8)

approve <- c(d_A, d_B, d_C, d_D, d_E)
industry <- as.factor(c(rep('A', 5),  rep('B', 10), rep('C', 20),
              rep('D', 30), rep('E', 40) ))

mod_mat <- model.matrix(lm(approve ~ industry))

d_list <- list(n = nrow(mod_mat),
               p = ncol(mod_mat),
               X = mod_mat, 
               approve = approve)

################################################################################
#####                       Bayesian Esimtate                              #####
################################################################################

mod<-stan_model(file="PartialPool.stan")

stan_res <- sampling(object = mod,
                     seed=11,data = d_list,
                     pars = c("p_approve"),
                     chains=1, iter=20000, warmup=10000)

plot(stan_res, plotfun='trace') ## convergence check

p_approve_sum <- summary(stan_res)

partial_pool_p <- data.frame(p_approve_sum$summary[-6,c('mean','2.5%','97.5%')])
partial_pool_p$industry = c('A','B','C','D','E')



################################################################################
#####                       Frequentist Estimates                          #####
################################################################################

### no pooling
no_pool_res <- glm(approve ~ industry, family = binomial(link = 'logit'))

no_pooled_p <- predict(no_pool_res, 
                    newdata = data.frame(industry=as.factor(c('A','B','C','D','E') )), 
                    type = 'response')

### complete pooling
pool_res <- glm(approve ~ 1, family = binomial(link = 'logit'))

pooled_p <- exp(pool_res$coefficients)/(1 + exp(pool_res$coefficients))

