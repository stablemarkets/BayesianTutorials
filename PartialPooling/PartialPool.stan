data {
  int n; // Number of observations
  int p;
  int approve[n];
  matrix[n, p] X;
}

parameters {
  vector[p] beta;
}

transformed parameters {
  vector[n] eta;
  eta = X * beta;
}

model {
  // specify priors
  beta ~ normal(0,3);

  // likelihood 
  approve  ~ bernoulli_logit(eta);
}

generated quantities{
  vector[p] p_approve;
  
  p_approve[1] = inv_logit(beta[1]);
  p_approve[2] = inv_logit(beta[1] + beta[2]);
  p_approve[3] = inv_logit(beta[1] + beta[3]);
  p_approve[4] = inv_logit(beta[1] + beta[4]);
  p_approve[5] = inv_logit(beta[1] + beta[5]);
}
