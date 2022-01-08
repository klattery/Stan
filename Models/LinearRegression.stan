// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower=0> P; 
  matrix[N, P] ind;
  vector[N] dep;
  vector[P] beta_prior;
  real<lower = 0> beta_prior_sigma;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[P] beta;
  real<lower = 0> beta_sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  vector[N] pred_ind = ind * beta;
  real intercept = mean(dep - pred_ind); // add mean of difference
  vector[N] pred = pred_ind + intercept;
  real ols_sig = sqrt(mean(square(dep - pred)));
  dep ~ normal(pred, ols_sig); // standard OLS regression
  beta ~ normal(0, beta_sigma); // ridge shrinkage
  beta_sigma ~ normal(beta_prior,beta_prior_sigma);
}

