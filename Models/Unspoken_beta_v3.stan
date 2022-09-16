// Kevin Lattery Nov 2021
// Conjoint Model in Stan
// con_sign has hard constraints for sign of coded parameter (pos/neg)
// paircon_matrix has pairwise constraints beyond con_sign

functions{
  matrix logistic_hinge(matrix x, matrix delta) {
    return delta .* log1p_exp(x ./ delta);
    // https://statmodeling.stat.columbia.edu/2017/05/19/continuous-hinge-function-bayesian-modeling/
  }
  
  real MNL_LL_par(array[] int array_slice,
                  int a_beg, int a_end, // Stan determines values of these for specific core
                  matrix beta_ind,  matrix X,  vector dep, array[] real wts,
                  array[] int start, array[] int end,
                  array[] int task_individual,
                  array[] int task_type,
                  vector task_timez_flip, vector beta_scale, 
                  vector conrate_conv, vector conrate_attr
  ) {
    real ll = 0;
    for (t in a_beg:a_end){
     if (task_type[t] == 1){ // Attraction
         real wt_task = wts[t] * 4/3 * (.25 + inv_logit(task_timez_flip[t] * conrate_attr[task_individual[t]]));
         ll+= dot_product(dep[start[t]:end[t]] * wt_task,
              log_softmax(X[start[t]:end[t]] * 
                col(beta_ind,task_individual[t]))
              );

     } else {
       if (task_type[t] == 2){ // Conversion also applies scale factor to utilities from attraction
         real wt_task =  wts[t] * 4/3 * (.25 + inv_logit(task_timez_flip[t] * conrate_conv[task_individual[t]]));
         ll+= dot_product(dep[start[t]:end[t]] * wt_task,
              log_softmax(X[start[t]:end[t]] * 
                (beta_scale[task_individual[t]] *
                col(beta_ind,task_individual[t]))
              ));       
        }
     }
    } // end for
    return ll;
  }
  
  matrix make_chol(array[] real d,array[] real tri_val,array[,] int tri_pos){
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
  
   array[,] int get_pos(matrix sym, int tri_n){
    int K = cols(sym);
    int pos = 1;
    array[tri_n, 2] int result;
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
  array[3] int<lower = 0> sizes; // 1=# columns coded, 2 = # cat vars w/levels, 3=# rows in code_master 
  matrix[N, sizes[1]] ind_coded; // Coded data mapped to ind 
  array[N, sizes[2]]int<lower = 0> ind_levels;

  // coding for each attribute combined in code_master, and code_blocks that define attributes in that
  matrix[sizes[3], P] code_master;
  int n_atts; // rows in code blocks = # atts
  array[n_atts, 5] int<lower = 0> code_blocks;
  
  // Upper level priors
  cov_matrix[P] prior_cov;  // Prior covariance matrix, recommend Lenk adjustment
  int<lower = 0> df; // df over P, recommend 5
  vector[P] prior_alpha;
  real<lower = 0> a_sig;  // alpha ~ N(prior_alpha, a_sig)
 
  // Constraints
  int<lower = 0, upper =1> con_use; // 0 = ignore constraints, 1 = use_constraints
  vector[P] con_sign; // Sign of constraints -1, +1 or 0
  int paircon_rows; // # rows in pairs constraint matrix, 0 for none
  matrix[paircon_rows, P] paircon_matrix; // Pair constraints: beta * paircon_matrix > 0
  
  // Covariates
  int P_cov; // Number of Covariate parameters (coded)
  matrix[I,P_cov] i_cov; // Resp covariates for each individual I
  
  // Misc: Weights, blocking, scale mult, multi-threading split
  array[T] real wts; // array of weights for each task
  matrix[P,P] cov_block; // Specifies blocks of covariance items
  real<lower = 0> prior_cov_scale;  // Multiply prior_cov by this:  Typically 1
  int<lower = 0> splitsize; // grainsize for parallel processing 
  
  // Ragged array matching, For each task in 1:T, Specify:
  array[T] int<lower = 1, upper = I> task_individual;
  array[T] int<lower = 1, upper = N> start;
  array[T] int<lower = 1, upper = N> end;
  
  // New Data for Unspoken not in Standard HB 
  vector[2] scale_bounds; // lb, ub for scale parameter attraction to conversion
  array[T] int task_type; // 1 = Attraction, 2 = Conversion
  vector[T] task_timez; // standardized time for the task within respondent x task_type x round
  real<lower = 0> timescale_con_hyper; // KL likes 2: conrate_conv ~ beta(1,hyper);
  real<lower = 0> timescale_att_hyper; // KL likes 2: conrate_attr ~ beta(1,hyper);
}

transformed data{
  vector[T] task_timez_flip = -1 * task_timez; // faster times will be bigger numbers
  
  matrix[N, P] ind; // ind_coded and ind_levels will expand to ind (below)
  matrix[P_cov, I] i_covt = i_cov'; // transpose of i_cov is what we use
  matrix[P,P] L = cholesky_decompose(prior_cov_scale * prior_cov/(P + df)); // For Wishart
  array[P] real df_chi;
  int tri_n = tri_sum(cov_block); // # of lower tri elements 
  array[tri_n,2]int tri_pos= get_pos(cov_block, tri_n); // position of elements

  int con_n = vec_not0(con_sign);
  array[con_n] int con_p;               // Array of which parameters are sign constrained
  matrix[con_n, I] con_delta;     // Con function scale for parameter and respondent
  int paircon_use = 0;
  
  array[T] int array_slice;  // Parallel threads slice across tasks 1:T  
  int count = 1;
  int lev_col = 1;
  for (i in 1:T){array_slice[i] = i;}
  for (i in 1:P){
    df_chi[i] = P + df - i + 1;
    if (fabs(con_sign[i]) > .00000001){
      con_p[count] = i;
      con_delta[count,1:I] = rep_row_vector(con_sign[i], I);
      count += 1;
    }
  }
  if (con_use == 0) con_n = 0; // set number of ordinal constraints to if no constraints
  if ((paircon_rows > 0) && (con_use == 1)) paircon_use = 1; // use paircon with rows and flag set
  // Coding of attributes with levels
  for (i in 1:n_atts){
    if (code_blocks[i,5] > 0){ //already coded, just copy
      int col_beg = code_blocks[i, 1];
      int col_end = code_blocks[i, 2];
      int col_coded = code_blocks[i, 5]; // beginning of coded column
      ind[:,col_beg:col_end] = ind_coded[:,col_coded:(col_coded + col_end - col_beg)];
    } else {
      int col_beg = code_blocks[i,1];
      int col_end = code_blocks[i,2];
      int row_beg = code_blocks[i,3];
      int row_end = code_blocks[i,4];
      int nrows = row_end - row_beg +1;
      int ncols = col_end - col_beg + 1;
      matrix[nrows, ncols] code_att =  code_master[row_beg:row_end, col_beg:col_end]; // code for this attribute
      row_vector[ncols] zeroes = rep_row_vector(0, ncols);
      for (n in 1:N){
        int lev = ind_levels[n,lev_col];
        if (lev >= 1) ind[n,col_beg:col_end] = code_att[lev,];
        if (lev == 0) ind[n,col_beg:col_end] = zeroes;
      } // end 1:N loop coding
      lev_col = lev_col + 1;
    } // end if  
  } // end for (att coding)  
} // end transformed data


parameters {
  vector[P] alpha; // upper MVN: mean vector of ind betas
  array[P] real<lower=0> bart_c ; // Bartlett diag^2
  array[tri_n]real bart_z; // Bartlett lower tri
  matrix[P, I] z; // individual z scores N(0,1)
  matrix[P, P_cov] i_cov_load; // loadings for i_cov (resp covariates)
  
  real<lower = scale_bounds[1], upper = scale_bounds[2]> beta_scale_mu; //mean
  real<lower = 0, upper = (scale_bounds[2] - scale_bounds[1])> beta_scale_sigma; // sd of scale across respondents
  vector<lower = scale_bounds[1], upper = scale_bounds[2]>[I] beta_scale;  // respondent specific scale factors
  vector<lower = 0, upper = 1>[I] timescale_con; // conversion: how informative time is -- conversion scale from timez
  vector<lower = 0, upper = 1>[I] timescale_att; // attraction: how informative time is -- conversion scale from timez
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
    if (con_n > 0){
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
  if (paircon_use == 1) target += -sum(log1p_exp(-100 * (paircon_matrix * beta_ind))); // penalty for soft constraints
  beta_scale ~ normal(beta_scale_mu,beta_scale_sigma);
  beta_scale_mu ~ normal(1, 2);
  beta_scale_sigma ~ normal(1, 2);
  timescale_con ~ beta(1,timescale_con_hyper); // biasing toward 0 = flat weights
  timescale_att ~ beta(1,timescale_att_hyper);
  target += reduce_sum(MNL_LL_par, array_slice, splitsize, 
                       beta_ind, ind, dep, wts, start, end, task_individual,
                       task_type, task_timez_flip, beta_scale, 
                       2 * timescale_con, 2 * timescale_att); // extend timescale to [0,2]
} // End Model
