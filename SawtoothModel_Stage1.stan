// Standard Conjoint Model
// Kevin Lattery Oct 2019
// Similar to Jim Savage code: https://gist.github.com/khakieconomics/0333a54dff4fabf6204080ca5bf92cb6

data {
  // Sizing Constants, Number of:
  int N; // Rows in data file
  int I; // Individuals
  int T; // Total unique tasks (individuals x mean tasks per individual)  
  int P; // Parameters = ncol(coded independent data)
  
  // Observed Data (Ind and Dep Variables)
  matrix[N, P] X; // Independent Data (coded)  
  vector<lower = 0, upper = 1>[N] choice; // Dep variable

  // Ragged array matching, For each task in 1:T, Specify:
  int<lower = 1, upper = I> task_individual[T]; // which Individual in 1:I does task
  int<lower = 1, upper = N> start[T]; // which row of data in [1:N] task begins
  int<lower = 1, upper = N> end[T]; // which row of data in [1:N] task ends
  
  cov_matrix[P] PriorCov;  // Prior covariance matrix
  // Recommend cov matrix with Lenk adjustment for # of levels
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
  cov_beta ~ inv_wishart(P + 5, PriorCov); // Others are possible
  // cov_beta ~ wishart(P + 2, PriorCov/(P + 2));
 
// log probabilities of each choice in the dataset
 {
   vector[N] loglike; // create a temporary holding vector - used to line 1 
   for(t in 1:T) {
    loglike[start[t]:end[t]] =
    log_softmax(X[start[t]:end[t]]*beta_ind[task_individual[t]]);       
     }
   target += dot_product(loglike,choice);
 }
}

