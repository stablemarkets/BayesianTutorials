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

n<-50 # number of observations
# simulate model matrix
x<-cbind(1, rnorm(n, 0, 1), rnorm(n, 5,10),rnorm(n, 100,10))

# true beta coefficients
tb<-c(1000, 50, -50, 10)

# true phi
tphi<-10000
I<-diag(1,n,n) # identity matrix used for covariance matrix

# simulate outcome for regression 
y<-t(rmvnorm(1, x%*%tb, tphi*I))

# simulate many outcomes...used later for asymptotic evaluations
y_list<-replicate(1000, t(rmvnorm(1, x%*%tb, tphi*I)),simplify = FALSE)

################################################################################
###### 1 - Run Blocked Gibbs Sampler
################################################################################

# function for blocked gibbs sampler
block_gibbs<-function(y, x, iter, burnin, trim){
  # initialize gibbs
  xprimex_inv<-solve(t(x)%*%x) # calculate once for repeated use in sampler
  phi<-numeric(iter) # shell for phi
  b<-matrix(nrow=iter, ncol = 4) # shell for betas
  phi[1]<-6 # initial phi value to start sampler
  
  # phi hyperparameters
  a<-.5
  g<-10000
  
  # gibbs sampling
  for(i in 2:iter ){
    b[i,]<-rmvnorm(n = 1, 
                   mean = ((xprimex_inv%*%t(x))%*%y), 
                   sigma = phi[i-1]*xprimex_inv )
    
    phi[i]<-rinvgamma(n = 1, 
                      shape = (n/2 + a), 
                      rate = .5*( t((y - x%*%t(t(b[i,])) ))%*%(y - x%*%t(t(b[i,])) ) ) + g)
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
  
  return(joint_post_long)
}

# run gibbs sampler with specified parameters
post_dist<-block_gibbs(y = y, x = x, iter = 500000, burnin = 100000, trim = 50)

################################################################################
###### 2 - Summarize and Visualize Posterior Distributions 
################################################################################

# calculate posterior summary statistics (stats not used in rest of code)
post_sum_stats<-post_dist %>%
  group_by(param) %>%
  summarise(median=median(draw),
            lwr=quantile(draw,.025),
            upr=quantile(draw,.975)) %>%
  mutate(true_vals=c(tb,tphi))

# merge on summary statistics
post_dist <- post_dist %>%
  left_join(post_sum_stats, by='param')

# plot MCMC Chains
ggplot(post_dist,aes(x=iter,y=draw)) +
  geom_line() +
  geom_hline(aes(yintercept=true_vals, col='red'), show.legend=FALSE)+
  facet_grid(param ~ .,scale='free_y',switch = 'y') +
  theme_bw() + 
  xlab('Gibbs Sample Iteration') + ylab('MCMC Chains') + 
  ggtitle('Gibbs Sampler MCMC Chains by Parameter')

# plot Posterior Distributions
ggplot(post_dist,aes(x=draw)) +
  geom_histogram(aes(x=draw),bins=50) +
  geom_vline(aes(xintercept = true_vals,col='red'), show.legend = FALSE) +
  facet_grid(. ~ param, scale='free_x',switch = 'y') +
  theme_bw() + 
  xlab('Posterior Distributions') + ylab('Count') + 
  ggtitle('Posterior Distributions of Parameters (true values in red)')

################################################################################
###### 3 - Assess Bias and Coverage
################################################################################

# run the estimation 1000 times to get 1000 posterior medians and CIs
bayes_res<-lapply(y_list, block_gibbs, x=x, iter=1000, burnin=100, trim=1)

calc_sumstats<-function(post_dist){
  post_sum_stats<-post_dist %>%
    group_by(param) %>%
    summarise(post_median=median(draw),
              lwr=quantile(draw,.025),
              upr=quantile(draw,.975)) %>%
    mutate(true_vals=c(tb,tphi))
  return(post_sum_stats)
}


all_sum_stats<-lapply(bayes_res, calc_sumstats)
all_sum_stats_stack<-bind_rows(all_sum_stats) %>%
  arrange(param) %>%
  rename(est=post_median)

eval_sum <- all_sum_stats_stack %>%
  mutate(covered=ifelse(true_vals<upr & true_vals>lwr,1,0)) %>%
  group_by(param) %>%
  summarise(est_var=var(est),
            est_mean=mean(est),
            bias=mean(est-true_vals),
            true_val=mean(true_vals),
            coverage=mean(covered)) %>%
  mutate(perc_bias=(bias/true_val)*100)

# format table
eval_sum<-eval_sum[,c(1,5,3,2,4,7,6)]
names(eval_sum)<-c('Parameter','True Value','Estimator Mean',
                   'Estimator Variance','Bias','Percent Bias (of truth)',
                   'Coverage of 95% CI')

print.xtable(xtable(eval_sum,
                    caption = 'Estimator Evaluation'),
             caption.placement = 'top', 
             include.rownames = F)
