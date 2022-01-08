// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower=0> P; 
  matrix[N, P] ind;
  vector[N] dep;
  vector[P] beta_mu;
  real<lower = 0> beta_sigma_hyper;
}

transformed data{
  vector [N] dep_z = (dep - mean(dep))/sd(dep); // standardize dep
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
  vector[N] pred_raw = ind * beta;
  vector[N] pred_z = (pred_raw - mean(pred_raw))/sd(pred_raw); // standardized 
  real ols_sig = sqrt(mean(square(dep_z - pred_z)));
  dep_z ~ normal(pred_z, ols_sig); // standard OLS regression
  beta ~ normal(beta_mu, beta_sigma); // ridge shrinkage
  beta_sigma ~ normal(0,beta_sigma_hyper);
}

