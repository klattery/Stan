// Standard Conjoint Model
// Kevin Lattery Oct 2019

data {
  // Sizing Constants, Number of:
  int N; // Rows in data file
  int P; // Parameters = ncol(coded independent data)
  int T; // Total unique tasks (individuals x mean tasks per individual)  
  int I; // Individuals

  // Observed Data (Ind and Dep Variables)
  matrix[N, P] X; // Independent Data (coded)  
  vector<lower = 0, upper = 1>[N] choice; // Dep variable

  cov_matrix[P] PriorCov;  // Prior covariance matrix - rec Lenk cov matrix
  int<lower = 0> df; // df over P, recommend 5
  
  // For each task in 1:T, Specify:
  int<lower = 1, upper = I> task_individual[T]; // which Individual in 1:I does task
  int<lower = 1, upper = N> start[T]; // which row of data in [1:N] task begins
  int<lower = 1, upper = N> end[T]; // which row of data in [1:N] task ends
}

transformed data{
  int df_tot = P + df; // Total degrees of freedom
  matrix[P, P] IW_S = PriorCov * (df_tot - 1); //Scale matrix for Inv Wishart
}

parameters {
  vector[P] alpha; // hypermeans of the part-worths
  cov_matrix[P] cov_beta; // cov of betas 
  vector[P] beta_ind[I]; // respondent level betas
}

model {
// priors on the parameters
  beta_ind ~ multi_normal(alpha, cov_beta);
  alpha ~ normal(0, 10);
  cov_beta ~ inv_wishart(df_tot, IW_S); 

// MNL log of predicted probabilities * choice
 {
   vector[N] logprob; // create a temporary holding vector 
   for(t in 1:T) {
    logprob[start[t]:end[t]] =
    log_softmax(X[start[t]:end[t]]*beta_ind[task_individual[t]]);       
     }
   target += dot_product(logprob,choice);
 }
}

