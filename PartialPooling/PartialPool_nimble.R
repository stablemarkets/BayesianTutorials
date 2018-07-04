library(nimble)
setwd("/Users/aoganisi/Dropbox/Stable Markets/BayesianTutorials/PartialPooling")

################################################################################
#####                       Simulate Data                                  #####
################################################################################
set.seed(10)

d_A<-rbinom(n = 10,  size = 1, prob =  .2)
d_B<-rbinom(n = 100, size = 1, prob = .3)
d_C<-rbinom(n = 10, size = 1, prob = .45)
d_D<-rbinom(n = 100, size = 1, prob = .7)
d_E<-rbinom(n = 10, size = 1, prob = .8)

approve <- c(d_A, d_B, d_C, d_D, d_E)
industry <- as.factor(c(rep('A', 10),  rep('B', 100), rep('C', 10),
              rep('D', 100), rep('E', 10) ))

mod_mat <- model.matrix(lm(approve ~ industry))

d_list <- list(X = mod_mat, 
               approve = approve)

p <- ncol(mod_mat)
n <- nrow(mod_mat)

################################################################################
#####                       Bayesian Esimtate                              #####
################################################################################

code <- nimbleCode({
  
  for(i in 1:p){
    beta[i] ~ dnorm(0, 3)
  }
  
  
  logit(eta[1:n]) <- X[1:n,1:p] %*% beta[1:p]
  
  
  for(i in 1:n) {
    approve[i] ~ dbern(prob = eta[i] )  
  }
  
  p_approve[1] <- expit(beta[1])
  p_approve[2] <- expit(beta[1] + beta[2])
  p_approve[3] <- expit(beta[1] + beta[3])
  p_approve[4] <- expit(beta[1] + beta[4])
  p_approve[5] <- expit(beta[1] + beta[5])
  
})

merge_model <- nimbleModel(code=code, 
                           constants=list(p=p, n=n),
                           inits = list(beta=c(0,0,0,0,0)),
                           data=d_list)

spec <- configureMCMC(merge_model)
spec$addSampler(type = 'RW_block', target ='beta',
                control = list(targetNodes='beta',
                               adaptive = TRUE ))
spec$monitors <- c('p_approve')



mcmc <- buildMCMC(spec)

compiled_model <- compileNimble(merge_model)
compiled_mcmc <- compileNimble(mcmc, project = merge_model)

compiled_mcmc$run(10000)

samples <- as.matrix(compiled_mcmc$mvSamples)
summary(samples)

partial_pool_p <- colMeans(samples[5000:10000,])

################################################################################
#####                       Frequentist Estimates                          #####
################################################################################

### no pooling
no_pool_res <- glm(approve ~ industry, family = binomial(link = 'logit'))

# compute probabilities of merger for each industry separately.
no_pooled_p <- predict(no_pool_res, 
                    newdata = data.frame(industry=as.factor(c('A','B','C','D','E') )), 
                    type = 'response')

### complete pooling
pool_res <- glm(approve ~ 1, family = binomial(link = 'logit'))

# compute probability of merger across all idustries, pooled.
pooled_p <- exp(pool_res$coefficients)/(1 + exp(pool_res$coefficients))


################################################################################
#####                       Visualize Results                              #####
################################################################################


plot(partial_pool_p, pch=20, col='red', ylim=c(0,1), axes=F, 
     xlab='Industry', ylab='Probability of Merger')
axis(1, at = 1:5, labels = paste0(unique(industry), " (n =",c(10,100,10,100,10),')' ) )
axis(2, at = seq(0,1,.2), labels= seq(0,1,.2) )

points(1:5, no_pooled_p, pch=20, col='black')
abline(h=pooled_p, lty=2)

legend('bottomright', 
       legend = c('Pooled Estimate','Stratified Estimates', 'Bayesian Estimate'),
       lty = c(2,NA,NA), col=c('black','black','red'), pch=c(NA, 20,20), bty='n')


