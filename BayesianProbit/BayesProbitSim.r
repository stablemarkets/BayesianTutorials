library(LaplacesDemon)
library(truncnorm)
library(mvtnorm)

set.seed(1)

n <- 2000

d_list <- replicate(n = 100, expr = sim_dat(n=n), simplify = F)

d <- d_list[[1]]


gibbs_iter <- 200000
burnin <- 10000
trimm <- 1
p <- 3

beta <- matrix(NA, nrow = p, ncol = gibbs_iter)
xm <- cbind(1, d$x1, d$x2)

z <- rnorm(n = n, mean = 0, sd = 1)

for(i in 1:gibbs_iter){
  
  beta[, i] <- rcond_post_beta(nparams = p, xm=xm, y = z, psi = 1, 
                                     beta_prior_mean = c(0,0,0), 
                                     beta_prior_var = c(10,10,10) )
  
  z[d$y==0] <- rtruncnorm(n = sum(d$y==0),a = -Inf, b = 0, mean = xm[d$y==0,] %*% beta[,i], sd = 1 )
  z[d$y==1] <- rtruncnorm(n = sum(d$y==1),a = 0, b = Inf, mean = xm[d$y==1,] %*% beta[,i], sd = 1 )
  
  
}

freqres <- summary(glm(data=d, formula = y ~ x1 + x2, family = binomial(link = 'probit')))
freqres$coefficients[,'Estimate']

post_draw_ind <- seq(burnin, gibbs_iter, trimm)
par(mfrow=c(3,1))
plot(beta[1,post_draw_ind], type='l')
abline(h=freqres$coefficients[1,'Estimate'], col='red', lwd=2)

plot(beta[2,post_draw_ind], type='l')
abline(h=freqres$coefficients[2,'Estimate'], col='red', lwd=2)

plot(beta[3,post_draw_ind], type='l')
abline(h=freqres$coefficients[3,'Estimate'], col='red', lwd=2)

