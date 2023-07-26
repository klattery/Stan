// Kevin Lattery April 2023
// Conjoint Model in Stan
// con_sign has hard constraints for sign of coded parameter (pos/neg)
// Speedup: More calcs inside multi-threading (Kees' idea)

functions{
  real MNL_LL_resp_slice(array[] int resp_slice,
                  int r_beg, int r_end, // Stan determines values for specific core
                  matrix ind,  vector dep, array[] real wts,
                  array[] int start, array[] int end,
                  array[] int task_individual,
                  array[,] int id_ranges, // id accounting
                  vector alpha, matrix i_cov_load, matrix i_covt, matrix cov_chol, matrix z, // beta stuff
                  int paircon_use, int con_n, array[] int con_p, matrix con_delta, matrix paircon_matrix,
                  matrix comb_probs,
                  vector ind_posts,vector ind_finite, vector ind_infinite,
                  vector b_finite, vector b_infinite_delta, vector none_posts,
                  vector dep_posts
  ) {
    real ll = 0;
    int n_resp = r_end - r_beg + 1;  // Number of respondents in slice

    // Estimate betas for slice of resondents
    matrix[rows(z),n_resp] beta_ind = rep_matrix(alpha, n_resp) + (i_cov_load * i_covt[,r_beg:r_end])+ cov_chol * z[,r_beg:r_end]; // Unconstrained
    if (con_n > 0){
      beta_ind[con_p,] = con_delta[,r_beg:r_end] .* log1p_exp(beta_ind[con_p,] ./ con_delta[,r_beg:r_end]);
    }
    if (paircon_use == 1) ll += -sum(log1p_exp(-100 * (paircon_matrix * beta_ind))); // soft constraints penalty
    
    // Loop through respondents and tasks for LL
    for (resp in r_beg:r_end){
      // id_ranges cols 1-5: task_beg, task_end, row_beg, row_end, n_rows
      vector[id_ranges[resp, 5]] util_all_exp = exp(ind[id_ranges[resp,3]:id_ranges[resp,4],] *
                                        col(beta_ind, resp - r_beg + 1)); // utilities across all resp tasks
      int row_adjust = id_ranges[resp,3]  - 1;  // Adjust row numbers in sliced data
      for (t in id_ranges[resp, 1]:id_ranges[resp, 2]){
        int t_start = start[t];
        int t_end = end[t];
        int n_alts = t_end - t_start + 1;
        vector[n_alts] util_task_exp  = util_all_exp[(t_start - row_adjust):(t_end - row_adjust)];
        ll+= dot_product(
                dep[t_start:t_end],
                log(comb_probs[t_start:t_end, 1:n_alts] *
                    util_task_exp/sum(util_task_exp)
                )) * wts[t]; // end dot_product, ll
        real utask_posts  = dot_product(util_task_exp, ind_posts[t_start:t_end]);
        if (utask_posts > 0){
          real u_sum = sum(util_task_exp);
          real u_fin = log(u_sum - dot_product(util_task_exp, ind_finite[t_start:t_end])) -
                       log(u_sum);
          real u_inf = log(u_sum - dot_product(util_task_exp, ind_infinite[t_start:t_end])) -
                       log(u_sum);
          real prob_posts = inv_logit(
            log(utask_posts) +
            b_finite[resp] * u_fin +
            (b_finite[resp] + b_infinite_delta[resp]) * u_inf -
            none_posts[resp]
          );
          ll+= log(prob_posts) * dep_posts[t] +
               log(1-prob_posts) * (1-dep_posts[t]);
        }
      } // end t loop through tasks
    } // end resp loop   
    return ll;
  }
  
  int vec_not0(vector vec){ // Count # of non-0 elements in vector
    int count = 0;
    for (i in 1:num_elements(vec)){
      if (abs(vec[i]) > .00000001) count += 1;
    }
    return count;
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
  array[N, sizes[2]]int<lower = 0> ind_levels; // Data to be coded

  // coding for each attribute combined in code_master, and code_blocks that define attributes in that
  matrix[sizes[3], P] code_master;
  int n_atts; // rows in code blocks = # atts
  array[n_atts, 5] int<lower = 0> code_blocks;
  
  // Upper level priors
  cov_matrix[P] prior_cov;  // Prior covariance matrix, recommend Lenk adjustment
  int<lower = 2> df; // df over P, recommend 2  - 5
  vector[P] prior_alpha;
  real<lower = 0> a_sig;  // alpha ~ N(prior_alpha, a_sig)
 
  // Constraints
  int<lower = 0, upper =1> con_use; // 0 = ignore constraints, 1 = use_constraints
  vector[P] con_sign; // Sign of constraints -1, +1 or 0
  int paircon_rows; // # rows in pairs constraint matrix, 0 for none
  matrix[paircon_rows, P] paircon_matrix; // Pair constraints: beta * paircon_matrix > 0
  array[2] real<lower = .01> con_factors; // multiplicative, bound
  
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
  array[I, 5] int<lower = 1, upper = N> id_ranges;  // task range (min, max), row_range of each id
  
  matrix[N, 7] comb_probs; // Combine Probs of Post and None
  vector[N] ind_posts; // Binary Indicator if aternative is optional Post
  vector[N] ind_finite; // Binary Indicator if aternative is finite post level
  vector[N] ind_infinite; // Binary Indicator if aternative is infinite post level
  vector[T] dep_posts; // Prob = % Max posts
}

transformed data{
  matrix[N, P] ind; // ind_coded and ind_levels will expand to ind (below)
  matrix[P_cov, I] i_covt = i_cov'; // transpose of i_cov is what we use
  int df_tot = P + df;
 // matrix[P,P] L = cholesky_decompose(prior_cov/df_tot); // For Wishart
  matrix[P,P] L_IW = cholesky_decompose(prior_cov * (df -1) * prior_cov_scale); // For Inverse Wishart

  int con_n = vec_not0(con_sign);
  array[con_n] int con_p;               // Array of which parameters are sign constrained
  matrix[con_n, I] con_delta;     // Con function scale for parameter and respondent
  int paircon_use = 0;
  vector[P] a_LB = rep_vector(-20, P); // LB for alpha (initial)
  vector[P] a_UB = rep_vector(20, P); // UB for alpha (initial)
  vector[P] prior_alpha_mod = prior_alpha; // Will modify alpha for constraints

  array[I] int resp_slice;  // Parallel threads slice across respondents 1:T  
  int count = 1;
  int lev_col = 1;
  for (r in 1:I){resp_slice[r] = r;} // Slice by respondents
  for (i in 1:P){
    if (abs(con_sign[i]) > .00000001){ 
      con_p[count] = i; // array of parameters that are constrained
      con_delta[count,1:I] = rep_row_vector(con_sign[i] * con_factors[1], I); // con_delta controls hinge function
      count += 1;
    }
  }
  if (con_use == 0) con_n = 0; // set number of ordinal constraints to 0 if no constraints
  if (con_n > 0){
    for (i in 1:P){
      if (con_sign[i] > 0){
        a_LB[i] = -con_factors[2]; // Pos Con LB alpha before hinge
        prior_alpha_mod[i] = -con_factors[2]; // lower bound will be mode
      } 
      if (con_sign[i] < 0){
        a_UB[i] = con_factors[2]; // Neg Con UB alpha before hinge
        prior_alpha_mod[i] = con_factors[2]; // upper bound will be mode
      } 
    }   
  }
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
  vector[P] alpha; // upper mean
  cholesky_factor_cov[P] cov_chol;
  matrix[P, I] z; // individual z scores N(0,1)
  matrix[P, P_cov] i_cov_load; // loadings for i_cov (resp covariates)
  vector<lower = 0>[I] b_finite;
  vector<lower = 0>[I] b_infinite_delta;
  vector[I] none_posts;
}

model {
  // priors on the parameters
  target+= normal_lupdf(alpha|prior_alpha_mod, a_sig);
  target += -sum(log1p_exp(-100 * (a_UB - alpha))); // Penalty if alpha > UB
  target += -sum(log1p_exp(-100 * (alpha - a_LB))); // Penalty if alpha < LB
  
  target+= inv_wishart_cholesky_lupdf(cov_chol|df_tot, L_IW);
  to_vector(z) ~ std_normal(); // respondent deviations from MVN(alpha, cov)
  if (P_cov > 0) to_vector(i_cov_load) ~ std_normal();
  b_finite ~ normal(1, 10);
  b_infinite_delta ~ normal(1, 10);
  none_posts ~ normal(0,10);
  target += reduce_sum(MNL_LL_resp_slice, resp_slice, splitsize,
                       ind, dep, wts, start, end, task_individual, id_ranges,
                       alpha, i_cov_load, i_covt, cov_chol, z,
                       paircon_use, con_n, con_p, con_delta, paircon_matrix,
                       comb_probs,
                       ind_posts,ind_finite, ind_infinite,
                       b_finite, b_infinite_delta, none_posts,
                       dep_posts);
} // End Model

generated quantities{
  matrix[P, I] beta_ind;
  {
  beta_ind = rep_matrix(alpha, I) + (i_cov_load * i_covt) + cov_chol * z; // Unconstrained
    if (con_n > 0){
      beta_ind[con_p,1:I] = con_delta .* log1p_exp(beta_ind[con_p,1:I] ./ con_delta);
    } 
  }
}

