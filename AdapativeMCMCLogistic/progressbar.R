
iter <- 1000000

pb<-txtProgressBar(min = 1, max = iter, style = 3)
t<-0
for(i in 1:iter){
  t<-t+i
  setTxtProgressBar(pb,value = i)
}



draws <- mvrnorm(n = 10, mu = rep(0,5), Sigma = diag(5))
cov(draws)
