// Sawtooth like - not efficient Conjoint Model in Stan
// Created for classroom use Kevin Lattery Sep 2022

data {
  // Sizing Constants, Number of:
  int N; // # Rows in data file
  int P; // # Parameters = ncol(coded independent data)
  int I; // # Respondents/Individuals
  int T; // Total unique tasks (individuals x mean tasks per individual)

  // Main data
  matrix[N, P] ind; // Independent data coded
  vector<lower = 0, upper = 1>[N] dep; // Dep variable

  // Upper level priors
  cov_matrix[P] prior_cov;  // Prior covariance matrix, recommend Lenk adjustment
  int<lower = 0> df; // df over P, recommend 5
  vector[P] prior_alpha;
  real<lower = 0> a_sig;  // alpha ~ N(prior_alpha, a_sig)
 
  // Misc: Weights, prior scale
  real<lower = 0> prior_cov_scale;  // Multiply prior_cov by this:  Typically 1
  array[T] real<lower = 0> wts; 
  
  // Ragged array matching, For each task in 1:T, Specify:
  array[T] int<lower = 1, upper = I> task_individual;
  array[T] int<lower = 1, upper = N> start;
  array[T] int<lower = 1, upper = N> end;
}

transformed data{
  int df_tot = P + df; 
  matrix[P,P] IW_S = prior_cov * df_tot; // For Wishart
}

parameters {
  vector[P] alpha; // hypermeans of the part-worths
  cov_matrix[P] cov_beta; // cov of betas 
  array[I] vector[P] beta_ind; // respondent level betas
}

model {
  // priors on the parameters
  beta_ind ~ multi_normal(alpha, cov_beta); 
  alpha ~ normal(prior_alpha, a_sig); // PriorAlpha can be vector of 0's or AggModel
  cov_beta ~ inv_wishart(df_tot, IW_S);
  for (t in 1:T){
    target += wts[t] * dot_product(
         log_softmax(ind[start[t]:end[t]] * beta_ind[task_individual[t]]),
                     dep[beta_ind[task_individual[t]]
         ); 
  }
} // End Model
