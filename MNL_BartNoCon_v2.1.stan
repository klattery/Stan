// Kevin Lattery June 2020
// Conjoint Model in Stan
// CovBlock defines blocks of parameters
// Wishart prior with Barlett Decomposition 
// Hinge Boundary Constraints
// Parallel threads using reduce_sum

functions{
  vector logistic_hinge(vector x, vector delta) {
    return delta .* log1p_exp(x ./ delta);
    // https://statmodeling.stat.columbia.edu/2017/05/19/continuous-hinge-function-bayesian-modeling/
  }
  
  real MNL_LL_par(int[] array_slice,
                  int a_beg, int a_end, // Stan determines values of these for specific core
                  matrix beta_ind,  matrix X,  vector dep,
                  int[] start, int[] end,
                  int[] task_individual
  ) {
    vector[end[a_end] - start[a_beg] +1] logprob; 
    int sub_adj = start[a_beg] -1; // logprob starts at 1, but first loop is start[a_beg]
    for (t in a_beg:a_end){
      logprob[(start[t]-sub_adj): (end[t]-sub_adj)]= 
        log_softmax(X[start[t]:end[t]] * 
                      col(beta_ind,task_individual[t]));
    }
    return dot_product(logprob,dep[start[a_beg]:end[a_end]]);
  }
  
  matrix make_chol(real[] d, real[] tri_val, int[ , ] tri_pos){
  // d is diagonal,
  // tri_val is vector of n elements for lower tri, tri_pos is nx2 array of positions
   int K = size(d);
   matrix[K, K] A = diag_matrix(to_vector(d)); 
   int count = 1;
   if (num_elements(tri_val) > 0){
     for (i in 1:size(tri_val)) {
       A[tri_pos[i,1], tri_pos[i,2]] = tri_val[count];
       count += 1;
     }
   }
   return(A);
  }

  int tri_sum(matrix sym){
  // get count of lower triangular elements  
   int K = cols(sym);
   int count = 0;
   for (i in 2:K){
     for (j in 1:(i-1)){
      if (fabs(sym[i,j]) > .00000001) count += 1;
     }
   }
   return(count);
  }

  int[,] get_pos(matrix sym, int tri_n){
    int K = cols(sym);
    int pos = 1;
    int result [tri_n, 2];
    for (i in 2:K){
      for (j in 1:(i-1)){
        if (fabs(sym[i,j]) > .00000001){
          result[pos,1] = i;
          result[pos,2] = j;
          pos += 1;
        }
      }
    }
    return(result);
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
  matrix[P,P] CovBlock; // Specifies blocks of items
  
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
  real df_chi[P];
  int tri_n = tri_sum(CovBlock); // # of lower tri elements 
  int tri_pos[tri_n,2] = get_pos(CovBlock, tri_n); // position of elements
  int array_slice[T]; // Parallel threads slice across tasks 1:T
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
  real bart_z [tri_n]; // Bartlett lower tri
  matrix[P, I] z; // individual z scores N(0,1)
}

transformed parameters {
  matrix[P, I] beta_ind;
  {
    matrix[P,P] cov_chol;
    if (tri_n == 0){
      // No off-diagonal
      cov_chol = L * diag_matrix(to_vector(sqrt(bart_c))); 
    } else {
     cov_chol = L * make_chol(sqrt(bart_c),bart_z,tri_pos); 
    }
    beta_ind = rep_matrix(alpha, I) + cov_chol * z;
  }
}

model {
  // priors on the parameters
  alpha ~ normal(0, 10); // Can also use PriorAlpha
  // for (i in 1:P) bart_c[i] ~ chi_square(P + df - i + 1);
  // bart_c ~ chi_square(df_chi);
  target+= chi_square_lpdf(bart_c|df_chi);
  if (tri_n > 0) bart_z ~ std_normal();
  to_vector(z) ~ std_normal(); // log probabilities of each choice in the dataset
  target += reduce_sum(MNL_LL_par, array_slice, splitsize, 
                       beta_ind, X, dep, start, end, task_individual);
} // End Model
