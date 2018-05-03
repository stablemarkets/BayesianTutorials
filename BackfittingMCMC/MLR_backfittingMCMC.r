################################################################################
###### 0 - Packages and Simulate Data
################################################################################
library(mvtnorm)
library(invgamma)
library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)
set.seed(200)

n<-100 # number of observation
# simulate model matrix
xvec <- rnorm(n, 0, 1.5)
x<-cbind(1, xvec, xvec^2,xvec^3)

# true beta coefficients
tb<-c(0, 50, -20, 10)

# true phi
tphi<-10000
I<-diag(1,n,n) # identity matrix used for covariance matrix

# simulate outcome for regression 
y<-t(rmvnorm(1, x%*%tb, tphi*I))
plot(xvec, y)
# simulate many outcomes...used later for asymptotic evaluations
y_list<-replicate(1000, t(rmvnorm(1, x%*%tb, tphi*I)),simplify = FALSE)

################################################################################
###### 1 - Run Blocked Gibbs Sampler
################################################################################

# function for blocked gibbs sampler
backfit_mcmc<-function(y, x, iter, burnin, trim){
  # initialize gibbs
  xprimex_inv<-solve(t(x)%*%x) # calculate once for repeated use in sampler
  phi<-numeric(iter) # shell for phi
  p <- ncol(x) # number of beta parameters
  b<-matrix(nrow=iter, ncol = p) # shell for betas
  pred_y <- matrix(nrow=length(y), ncol=iter)
  
  phi[1]<-6 # initial phi value to start sampler
  b[1,] <- rnorm(p) # random initial values for betas
  
  # phi hyperparameters
  a<-.5
  g<-10000
  
  # beta hyperparameters
  # mu_0 <- 0
  # phi_0 <- 1000
  
  # gibbs sampling
  for(i in 2:iter ){
    
    for(par in 1:p){
      # compute residuals after applying other parameters
      x_sub <- x[,-par]
      b_sub <- b[i, -par]
      b_sub_prev <- b[i-1, -par]
      b_sub[is.na(b_sub)] <- b_sub_prev[is.na(b_sub)]
      
      r <- y - x_sub %*% t(t(b_sub))

      # mean for this paramater
      x_j <- x[,par, drop=F]
      sum_xj_sq <- sum(x_j^2)
      sum_r_xj <- sum(as.vector(r)*as.vector(x_j))

      b[i, par]<-rnorm(n = 1, mean = sum_r_xj/sum_xj_sq , sd = sqrt( phi[i-1]/sum_xj_sq  ) )
    }

    
    phi[i]<-rinvgamma(n = 1, 
                      shape = (n/2 + a), 
                      rate = .5*( t((y - x%*%t(t(b[i,])) ))%*%(y - x%*%t(t(b[i,])) ) ) + g)
    
    pred_y[,i] <- rmvnorm(1, x%*% t(t(b[i,])), phi[i]*I ) # posterior predictive draw.
  }
  
  # apply burnin and trimming  
  keep_draws<-seq(burnin,iter,trim)
  phi<-phi[keep_draws]
  b<-b[keep_draws,]
  
  # format and output
  joint_post<-data.frame(b=b,phi=phi)
  colnames(joint_post)[1:(ncol(x))]<-paste0('B',0:(ncol(x)-1) )
  
  joint_post_long<-gather(joint_post,keep_draws) %>%
    rename(param=keep_draws, draw=value) %>%
    mutate(iter=rep(keep_draws,ncol(joint_post)))
  
  return(list(joint_post_long, pred_y))
}

# run gibbs sampler with specified parameters
post_dist<-backfit_mcmc(y = y, x = x, iter = 10000, burnin = 5000, trim = 1)


tt <- post_dist[[2]] # extract posterior predictive
png(filename = 'OriginalData.png')
plot(xvec, y, 
     xlab = 'x', ylab='y', 
     main='Simulated (x,y) Data', 
     ylim=c(-2500,1500), pch=20)
dev.off()

png(filename = 'BayesianResults.png')
plot(xvec, y, 
     xlab = 'x', ylab='y', 
     main='Posterior Mean Function fit using Backfitting MCMC', 
     ylim=c(-2500,1500), pch=20)
for(i in 5000:10000) lines(sort(xvec), sort(tt[,i]), col='gray')
lines(sort(xvec), sort(rowMeans(tt[,9000:10000]) ), col='red', lwd=3)
points(xvec, y, pch=20 )
legend('bottomright', 
       legend = c('Posterior Predictive Mean', '5000 posterior predictive draws','Data'),
       col = c('red','gray','black'), pch = c(NA, 22, 20), bty='n', 
       lty=c(1, NA, NA), lwd=c(1,5, NA)  )
dev.off()