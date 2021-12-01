// Kevin Lattery Nov 2021
// Conjoint Model in Stan
// con_sign has hard constraints for sign of coded parameter (pos/neg)
// paircon_matrix has pairwise constraints beyond con_sign

functions{
  matrix logistic_hinge(matrix x, matrix delta) {
    return delta .* log1p_exp(x ./ delta);
    // https://statmodeling.stat.columbia.edu/2017/05/19/continuous-hinge-function-bayesian-modeling/
  }
  
  real MNL_LL_par(int[] array_slice,
                  int a_beg, int a_end, // Stan determines values of these for specific core
                  matrix beta_ind,  matrix X,  vector dep,
                  int[] start, int[] end,
                  int[] task_individual
  ) {
    real ll = 0; 
    for (t in a_beg:a_end){
        ll+= dot_product(dep[start[t]:end[t]],
             log_softmax(X[start[t]:end[t]] * col(beta_ind,task_individual[t])));  
      }
    return ll;
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
   return A;
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
   return count;
  }
  
  int vec_not0(vector vec){
    int count = 0;
    for (i in 1:num_elements(vec)){
      if (fabs(vec[i]) > .00000001) count += 1;
    }
    return count;
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
  int P; // Parameters = ncol(coded independent data)
  int I; // Individuals
  int T; // Total unique tasks (individuals x mean tasks per individual)

  // Main data
  vector<lower = 0, upper = 1>[N] dep; // Dep variable
  int<upper = P> P_comp; // # parameters in compacted data 
  matrix[N, P_comp] ind_comp; // Independent Data (compact) becomes ind[N, P] 
  
  // for unpacking ind_comp
  int code_matrix_rows; // rows in code matrix
  matrix[code_matrix_rows, P] code_matrix;
  int<lower = 0> code_blocks[P_comp, 4]; // 1 = r_beg, 2 = r_end, 3 = col_beg, 4 = col_end
  
  // Upper level priors
  cov_matrix[P] prior_cov;  // Prior covariance matrix, recommend Lenk adjustment
  int<lower = 0> df; // df over P, recommend 5
  vector[P] prior_alpha;
  real<lower = 0> a_sig;  // alpha ~ N(prior_alpha, a_sig)
 
  // Constraints
  vector[P] con_sign; // Sign of constraints -1, +1 or 0
  int paircon_rows; // # rows in pairs constraint matrix, 0 for none
  matrix[paircon_rows, P] paircon_matrix; // Pair constraints: beta * paircon_matrix > 0
  
  // Covariates
  int P_cov; // Number of Covariate parameters (coded)
  matrix[I,P_cov] i_cov; // Resp covariates for each individual I
  
  // Misc: Weights, blocking, scale mult, multi-threading split
  vector[N] wts; // weights of each row
  matrix[P,P] cov_block; // Specifies blocks of covariance items
  real<lower = 0> prior_cov_scale;  // Multiply prior_cov by this:  Typically 1
  int<lower = 0> splitsize; // grainsize for parallel processing 
  
  // Ragged array matching, For each task in 1:T, Specify:
  int<lower = 1, upper = I> task_individual[T]; // which i in 1:I does task
  int<lower = 1, upper = N> start[T]; // which row in [1:N] task begins
  int<lower = 1, upper = N> end[T];   // which row in [1:N] task ends
}

transformed data{
  matrix[N, P] ind; // Fill in the ind
  matrix[P_cov, I] i_covt = i_cov'; // transpose of i_cov is what we use
  vector[N] dep_wt = dep .* wts;
  matrix[P,P] L = cholesky_decompose(prior_cov_scale * prior_cov/(P + df)); // For Wishart
  real df_chi[P];
  int tri_n = tri_sum(cov_block); // # of lower tri elements 
  int tri_pos[tri_n,2] = get_pos(cov_block, tri_n); // position of elements

  int con_n = vec_not0(con_sign); // Number of parameters constrained
  int con_p[con_n];               // Array of which parameters are constrained
  matrix[con_n, I] con_delta;     // Con function scale for parameter and respondent
  
  int array_slice[T]; // Parallel threads slice across tasks 1:T  
  int count = 1;
  int ind_beg = 1;
  int ind_end = 1;
  for (i in 1:P){
    df_chi[i] = P + df - i + 1;
    if (fabs(con_sign[i]) > .00000001){
      con_p[count] = i;
      con_delta[count,1:I] = rep_row_vector(con_sign[i], I);
      count += 1;
    }
  }
  for (i in 1:T){
    array_slice[i] = i;
  }
  
  for (i in 1:P_comp){ 
    if (code_blocks[i,1] == 0){ //no coding, just copy
      ind_end = ind_beg;
      ind[:,ind_beg:ind_end] = ind_comp[:,i:i];
    } else {
      int row_beg = code_blocks[i,1];
      int row_end = code_blocks[i,2];
      int col_beg = code_blocks[i,3];
      int col_end = code_blocks[i,4];
      int nrows = row_end - row_beg +1;
      int ncols = col_end - col_beg + 1;
      matrix[nrows, ncols] code_att =  code_matrix[row_beg:row_end, col_beg:col_end]; // code for this attribute
      ind_end = ind_beg + ncols -1;
      for (n in 1:N){
        real x_raw = ind_comp[n,i]; // level in the data
        int x_int = 0; // integer to convert xraw to
        while (x_int + .01 <  x_raw) x_int +=1;
        ind[n,ind_beg:ind_end] = code_att[x_int,];
      } // end loop coding
    } // end if  
    ind_beg = ind_end + 1; // Start at next column for ind output  
  } // end for p_comp 
} // end transformed data


parameters {
  vector[P] alpha; // upper MVN: mean vector of ind betas
  real<lower=0> bart_c [P]; // Bartlett diag^2
  real bart_z [tri_n]; // Bartlett lower tri
  matrix[P, I] z; // individual z scores N(0,1)
  matrix[P, P_cov] i_cov_load; // loadings for i_cov (resp covariates)
}

transformed parameters {
  matrix[P, I] beta_ind;
  {
    matrix[P,P] cov_chol;
    if (tri_n == 0){ // No off-diagonal
      cov_chol = diag_post_multiply(L, to_vector(sqrt(bart_c))); // L * diag_matrix()
    } else {
      cov_chol = L * make_chol(sqrt(bart_c),bart_z,tri_pos); 
    }
    beta_ind = rep_matrix(alpha, I) + (i_cov_load * i_covt) + cov_chol * z; // Unconstrained
    if (con_n >0){
      beta_ind[con_p,1:I] = con_delta .* log1p_exp(beta_ind[con_p,1:I] ./ con_delta);
    } 
  }
}

model {
  // priors on the parameters
  alpha ~ normal(prior_alpha, a_sig); // PriorAlpha can be vector of 0's or AggModel
  target+= chi_square_lpdf(bart_c|df_chi); // for (i in 1:P) bart_c[i] ~ chi_square(P + df - i + 1);
  if (tri_n > 0) bart_z ~ std_normal();
  to_vector(z) ~ std_normal(); // log probabilities of each choice in the dataset
  if (P_cov > 0) to_vector(i_cov_load) ~ std_normal();
  if (paircon_rows > 0) target += -sum(log1p_exp(-100 * (paircon_matrix * beta_ind))); // penalty for soft constraints
  target += reduce_sum(MNL_LL_par, array_slice, splitsize, 
                       beta_ind, ind, dep_wt, start, end, task_individual);
} // End Model
