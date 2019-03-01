#Author: Arman Oganisian

library(microbenchmark)
library(LaplacesDemon)
library(invgamma)
library(MASS)
library(profvis)
source("HelperFunctions.R")

set.seed(10)


################################################################################
### 1 - Run Benchmark
################################################################################
# hyper-parameters
lambda<-c(0,0,0,0)
phi<-10

n_vec <- seq(1000, 16000, 2000)
rel_time <- numeric(length = length(n_vec))

for(i in 1:length(n_vec)){
  
  N<-n_vec[i]
  d<-data.frame(age_group=sample(x = c(0,1,2), size = N, replace = T))
  d$age_1<-ifelse(d$age_group==1,1,0)
  d$age_2<-ifelse(d$age_group==2,1,0)
  
  d$trt<-rbinom(n = N, size = 1,prob = invlogit(0 + 2*d$age_1 + - 2*d$age_2))
  
  d$y<-rbinom(n = N, size = 1,
              prob = invlogit(-1 + .7*d$age_1 + 1.1*d$age_2 + 1.1*d$trt))
  
  X<-as.matrix(cbind(1,d[,2:4])) # model matrix
  Y<-matrix(d$y, ncol=1) # outcome vector
  
  bench<-microbenchmark(MH_vanilla = mh_vanilla(beta_0 = c(0,0,0,0), 
                                                phi = phi, lambda = lambda, X = X, Y = Y, 
                                                mh_trials=1000, jump_v=.2),
                        MH_cache = mh_cache(beta_0 = c(0,0,0,0), 
                                            phi = phi, lambda = lambda, X = X, Y = Y,
                                            mh_trials=1000, jump_v=.2),
                        times = 10)
  
  bench_sum <- summary(bench)
  vanilla_time <- bench_sum$mean[bench_sum$expr=='MH_vanilla']
  cache_time <- bench_sum$mean[bench_sum$expr=='MH_cache']
  rel_time[i] <- vanilla_time/cache_time
}

################################################################################
### 3 - Plot Results
################################################################################

par(mfrow=c(1,1))
plot(n_vec, rel_time, type ='l')
