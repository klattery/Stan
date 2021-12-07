// Standard Conjoint Model
// Kevin Lattery Oct 2019
// Similar to Jim Savage code: https://gist.github.com/khakieconomics/0333a54dff4fabf6204080ca5bf92cb6

functions{
  real log_ratio(matrix m_test, matrix m_base, real n_test, real n_base){
    real n = n_test + n_base;
    real p1 = n_test * (1/n);
    real p2 = n_base * (1/n);
    matrix[cols(m_test),cols(m_test)] m_mean = (p1 * m_test) + (p2 * m_base);
    return (n_test/2) * log_determinant(m_test) +
           (n_base/2) * log_determinant(m_base) - 
           (n/2) * log_determinant(m_mean);
    }
}

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
  corr_matrix[P] PriorCor; // Converted to Corr Matrix
  // Recommend cov matrix with Lenk adjustment for # of levels
}

parameters {
  vector[P] alpha; // hypermeans of the part-worths
  corr_matrix[P] cor_beta; // prior correlation
  vector<lower=0>[P] tau; // prior scale
  matrix[I, P] z; // individual z scores N(0,1)
}

transformed parameters {
  // beta_ind ~ MVN(alpha, cov_beta) iff beta_ind = alpha + z * sigma
  matrix[I, P] beta_ind = rep_matrix(alpha', I) +
    z * (cholesky_decompose(quad_form_diag(cor_beta, tau)));
}

model {
// priors on the parameters
  to_vector(z) ~ normal(0, 1); // z is ~ z scores
  alpha ~ normal(0, 10);
  tau ~ normal(1, 5);
  
  // cor_beta ~ lkj_corr(4); // No Prior as Input
  // cor_beta ~ inv_wishart(P + 2, PriorCor);
  target += log_ratio(cor_beta, PriorCor, N, P);
  
// log probabilities of each choice in the dataset
 {
   vector[N] loglike; // create a temporary holding vector - used to line 1 
   for(t in 1:T) {
    loglike[start[t]:end[t]] =
    log_softmax(X[start[t]:end[t]] * beta_ind[task_individual[t]]');       
     }
   target += dot_product(loglike,choice);
 }
}

