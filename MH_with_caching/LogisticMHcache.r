#Author: Arman Oganisian

library(microbenchmark)
library(LaplacesDemon)
library(MASS)
source("HelperFunctions.R")

set.seed(10)


################################################################################
### 1 - Run Benchmark
################################################################################
# hyper-parameters and true values
lambda<-c(0,0,0,0)
phi<-10
true_beta <- matrix(c(0,2,1,-2),ncol=1)
N<-50000

# simulate covariates 
X1 <- X2 <-X3 <- rnorm(N)
X <- model.matrix(~ X1 + X2 + X3)

# simulate outcome
Y <- rbinom(n = N, size = 1, prob = invlogit( X %*% true_beta  ) )


bench<-microbenchmark(
                      # Run Vanialla Metropolis
                      MH_vanilla = mh_vanilla(beta_0 = c(0,0,0,0), # initial value
                                              p=4, # number of parameters
                                              phi = phi, lambda = lambda, #hyperparameters
                                              X = X, Y = Y, # Data
                                              #iterations and proposal variance
                                              mh_trials=1000, jump_v=.2 ),
                      # Run Metropolis with Cache
                      MH_cache = mh_cache(beta_0 = c(0,0,0,0), 
                                          phi = phi, lambda = lambda, X = X, Y = Y,
                                          mh_trials=1000, jump_v=.2, p=4),
                      times = 10)
bench

