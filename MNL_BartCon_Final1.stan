// Conjoint Model in Stan
// Hinge Boundary Constraints
// Wishart prior with Barlett Decomposition 
// Kevin Lattery Apr 2020
// Parallel threads using reduce_sum

functions{
  vector logistic_hinge(vector x, vector delta) {
    return delta .* log1p_exp(x ./ delta);
    // https://statmodeling.stat.columbia.edu/2017/05/19/continuous-hinge-function-bayesian-modeling/
  }
  
  real MNL_LL(vector[] beta_ind, matrix X, vector dep,
              int N, int T, int[] start, int[] end, int[] task_individual){
    vector[N] logprob; // create a temporary holding vector  
    for(t in 1:T) {
      logprob[start[t]:end[t]] =
        log_softmax(X[start[t]:end[t]] * beta_ind[task_individual[t]]);     
    }
    return dot_product(logprob,dep);
  }
  
  real MNL_LL_par(int[] array_slice,
                  int a_beg, int a_end, // Stan determines values of these for specific core
                  vector[] beta_ind,  matrix X,  vector dep,
                  int[] start, int[] end,
                  int[] task_individual
  ) {
    vector[end[a_end] - start[a_beg] +1] logprob; 
    int sub_adj = start[a_beg] -1; // logprob starts at 1, but first loop is start[a_beg]
    for (t in a_beg:a_end){
      logprob[(start[t]-sub_adj): (end[t]-sub_adj)]= 
        log_softmax(X[start[t]:end[t]] * 
                      beta_ind[task_individual[t]]);
    }
    return dot_product(logprob,dep[start[a_beg]:end[a_end]]);
  }
  
  matrix make_lower(real[] d, real[] low_tri){
  // d is diagonal, low_tri is vector of dC2 elements for lower tri
  int K = size(d);
  int count;
  matrix[K, K] A;
  count = 1;
  for (j in 1:(K-1)) {
    for (i in (j+1):K) {
      A[i, j] = low_tri[count];
      count += 1;
    }
    for (i in 1:(j - 1)) A[i, j] = 0.0;
    A[j, j] = d[j];
  }
  for (i in 1:(K-1)) A[i, K] = 0;
  A[K, K] = d[K];
  return(A);
  }
} // End Functions

data {
  // Sizing Constants, Number of:
  int N; // Rows in data file
  int I; // Individuals
  int T; // Total unique tasks (individuals x mean tasks per individual)  
  int P; // Parameters = ncol(coded independent data)
  
  // Observed Data (Ind and Dep Variables)
  matrix[N, P] X; // Independent Data (coded)  
  vector<lower = 0, upper = 1>[N] choice; // Dep variable
  
  cov_matrix[P] PriorCov;  // Prior covariance matrix, recommend Lenk adjustment
  vector[P] PriorAlpha;
  vector[P] con_sign; // Sign of constraints -1, +1 or 0
  vector[N] weights; // weights of each row
  
  // Ragged array matching, For each task in 1:T, Specify:
  int<lower = 1, upper = I> task_individual[T]; // which i in 1:I does task
  int<lower = 1, upper = N> start[T]; // which row in [1:N] task begins
  int<lower = 1, upper = N> end[T];   // which row in [1:N] task ends
  int splitsize; // grainsize for parallel processing
  // matrix[P,P] CovBlock;
  int<lower = 0> df; // df over P, recommend 5
}

transformed data{
  vector[N] dep = choice .* weights;
  vector[P] bounded; // binary as to bounded or not
  vector[P] delta; // scale factor for bounding function
  matrix[P,P] L = cholesky_decompose(PriorCov/(P + df)); // For Wishart
  int array_slice[T]; // Parallel threads slice across tasks 1:T
  real df_chi[P];
  for (i in 1:T){
    array_slice[i] = i;
  }  
  for (i in 1:P){
    df_chi[i] = P + df - i + 1;
    if (con_sign[i] == 0){
      bounded[i] = 0; // Not bounded
      delta[i] = 1; // Any non-zero value
    } else {
      bounded[i] = 1; // Bounded
      delta[i] = con_sign[i];
    }
  }   
}

parameters {
  vector[P] alpha; // upper MVN: mean vector of ind betas
  real<lower=0> bart_c [P]; // Bartlett diag^2
  real bart_z [P * (P - 1)/2]; // Bartlett lower tri
  matrix[P, I] z; // individual z scores N(0,1)
}

transformed parameters {
  vector[P] beta_ind[I];
  {
    matrix[P,P] cov_chol = L * make_lower(sqrt(bart_c),bart_z); 
    matrix[P, I] beta_raw = rep_matrix(alpha, I) + cov_chol * z;
    for (i in 1:I){
      beta_ind[i] = ((1- bounded) .* beta_raw[,i]) + 
        (bounded .* logistic_hinge(beta_raw[,i], delta));
    }  
  }
}

model {
  // priors on the parameters
  alpha ~ normal(0, 10); // Can also use PriorAlpha
  // for (i in 1:P) bart_c[i] ~ chi_square(P + df - i + 1);
  // bart_c ~ chi_square(df_chi);
  target+= chi_square_lpdf(bart_c|df_chi);
  bart_z ~ std_normal();
  to_vector(z) ~ std_normal(); // log probabilities of each choice in the dataset
  target += reduce_sum(MNL_LL_par, array_slice, splitsize, 
                       beta_ind, X, dep, start, end, task_individual);
} // End Model
