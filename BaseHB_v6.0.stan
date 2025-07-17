// Kevin Lattery April 2023
// Conjoint Model in Stan
// con_sign has hard constraints for sign of coded parameter (pos/neg)
// Speedup: More calcs inside multi-threading (Kees' idea)
// Speedup by sorting chosen alts on top

functions{
  real MNL_LL_resp_slice(array[] int resp_slice,
                  int r_beg, int r_end, // Stan determines values for specific core
                  matrix ind_mult,  vector dep, array[] real wts, array[] int holdout,
                  array[] int start, array[] int end,
                  array[] int npos_m1, array[] int task_individual,
                  array[,] int id_ranges, // id accounting
                  matrix code_master_block, int n_atts_pick,
                  array[,] int beta_full_pick, array[] int cols_mult_pos, // for back coding
                  vector alpha, matrix i_cov_load, matrix i_covt, matrix cov_chol, matrix z, // beta stuff
                  int paircon_use, int con_n, array[] int con_p, matrix con_delta, matrix paircon_matrix,
                  array[] int con_p_pos, array[] int con_p_neg, array[] real con_factors,
                  int check_bad, array[] real prob_bad, array[] real ll0_resp,
                  array[] real scale_factor,array[] int task_scale_group
  ) {
    real ll = 0;
    int n_resp = r_end - r_beg + 1;  // Number of respondents in slice
    int P_full_p1 = rows(code_master_block) + 1; // # Back coded betas + 1
    // int check_bad = num_elements(prob_bad); // 0 or 1
    real log_good = 0;
    real log_bad = 0;
    if (check_bad == 1){
      log_bad = log(prob_bad[1]);
      log_good = log(1 - prob_bad[1]);
    }

    // Estimate betas for slice of resondents
    matrix[rows(z),n_resp] beta_ind;
    if (con_n > 0){
      matrix[rows(z),n_resp] beta_dev = (i_cov_load * i_covt[,r_beg:r_end])+ cov_chol * z[,r_beg:r_end]; // Deviations from mean
      if (num_elements(con_p_neg) > 0) ll+= -(100/con_factors[4]) * sum(log1p_exp(con_factors[3] * (con_factors[4] - beta_dev[con_p_neg,]))); // Penalty if con_sign deviation > dev_cap
      if (num_elements(con_p_pos) < 0) ll+= -(100/con_factors[4]) * sum(log1p_exp(con_factors[3] * (beta_dev[con_p_pos,] - con_factors[4]))); // Penalty if con_sign deviation < -dev_cap
      beta_ind = rep_matrix(alpha, n_resp) + beta_dev;
      beta_ind[con_p,] = con_delta[,r_beg:r_end] .* log1p_exp(beta_ind[con_p,] ./ con_delta[,r_beg:r_end]);
    } else {
      beta_ind = rep_matrix(alpha, n_resp) + (i_cov_load * i_covt[,r_beg:r_end])+ cov_chol * z[,r_beg:r_end]; // Unconstrained
    }
    if (paircon_use == 1) ll += -sum(log1p_exp(con_factors[5] * (paircon_matrix * beta_ind + con_factors[6]))); // soft constraints penalty
    
    // Loop through respondents and tasks for LL
    int n_cols_mult = num_elements(cols_mult_pos);
    for (resp in r_beg:r_end){
      vector[id_ranges[resp, 5]] util_all;
      // id_ranges cols 1-5: task_beg, task_end, row_beg, row_end, n_rows
      if (n_atts_pick == 0){
        util_all = ind_mult[id_ranges[resp,3]:id_ranges[resp,4],] *
                                        col(beta_ind, resp - r_beg + 1); // utilities across all resp tasks
      } else {
       util_all = rep_vector(0,id_ranges[resp, 5]); // Build utility since we have atts with levels to pick
       if (n_cols_mult > 0){
          util_all += ind_mult[id_ranges[resp,3]:id_ranges[resp,4],] *
                      col(beta_ind, resp - r_beg + 1)[cols_mult_pos]; // utilities across all resp tasks
       }
       vector[P_full_p1] beta_full = append_row(0, 
                 code_master_block * col(beta_ind, resp - r_beg + 1));
                 // full set of back coded utilities with 0 prefix
       for (att_pick in 1:n_atts_pick){
         util_all += beta_full[beta_full_pick[id_ranges[resp,3]:id_ranges[resp,4], att_pick]];
       }  
      }
      // LL for resp across tasks
      int row_adjust = id_ranges[resp,3]  - 1;  // Adjust row numbers in sliced data
      real ll_resp = 0; // holder to sum LL across tasks
      for (t in id_ranges[resp, 1]:id_ranges[resp, 2]){
        if (holdout[t] == 0) {
        int start_row = start[t] - row_adjust;
        if (task_scale_group[t]>1){
              util_all[start_row:(end[t] - row_adjust)]=util_all[start_row:(end[t] - row_adjust)]*scale_factor[task_scale_group[t]-1];
        }
        ll_resp+=  dot_product(
                  dep[start[t]:  (start[t] + npos_m1[t])],
                  util_all[start_row:(start_row + npos_m1[t])] -
                  log_sum_exp(util_all[start_row:(end[t] - row_adjust)])
                ) * wts[t];
        }
      } // end t loop through tasks
      if (check_bad == 0){
        ll+= ll_resp;
      } else {
        ll+=log_sum_exp(log_good + ll_resp, log_bad + ll0_resp[resp]);
      }
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
  int vec_pos(vector vec){ // Count # of non-0 elements in vector
    int count = 0;
    for (i in 1:num_elements(vec)){
      if (vec[i] > .00000001) count += 1;
    }
    return count;
  }
  int vec_neg(vector vec){ // Count # of non-0 elements in vector
    int count = 0;
    for (i in 1:num_elements(vec)){
      if (vec[i] < -.00000001) count += 1;
    }
    return count;
  }
// Functions fo Bartlett Decomposition
  matrix make_chol(array[] real d, array[] real tri_val, array[,] int tri_pos){
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
      if (abs(sym[i,j]) > .00000001) count += 1;
     }
   }
   return count;
  }
  
  array[,] int get_pos(matrix sym, int tri_n){
    int K = cols(sym);
    int pos = 1;
    array[tri_n, 2] int result;
    for (i in 2:K){
      for (j in 1:(i-1)){
        if (abs(sym[i,j]) > .00000001){
          result[pos,1] = i;
          result[pos,2] = j;
          pos += 1;
        }
      }
    }
    return(result);
  }

  array[,] int get_pos_bart(int P, int n_cols_bart, int tri_n){
    int pos = 1;
    array[tri_n, 2] int result;
    for (i in 2:P){
      for (j in 1:(i-1)){
        if (j <= n_cols_bart){
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
  // int P_fac;  // Factor reduced # of Paramters
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
  array[6] real con_factors; // multiplicative, bound, dev_mult, dev_cap, con_pair
  // real<upper = 0> con_pair_penalty;
  
  // Covariates
  int P_cov; // Number of Covariate parameters (coded)
  matrix[I,P_cov] i_cov; // Resp covariates for each individual I
  
  // Misc: Weights, blocking, scale mult, multi-threading split
  array[T] real wts; // array of weights for each task
  // matrix[P,P] cov_block; // Specifies blocks of covariance items
  real<lower = 0> prior_cov_scale;  // Multiply prior_cov by this:  Typically 1
  int<lower = 1> splitsize_n;  // 1 means we are estimating, more than 1 means testing
  array[splitsize_n] int<lower = 0> splitsize; // grainsize for parallel processing 
  
  // Ragged array matching, For each task in 1:T, Specify:
  array[T] int<lower = 1, upper = I> task_individual;
  array[T] int<lower = 1, upper = N> start;
  array[T] int<lower = 1, upper = N> end;
  array[T] int<lower = 0> npos_m1; // Number of choices > 0 each task minus 1
  array[I, 5] int<lower = 1, upper = N> id_ranges;  // task range (min, max), row_range of each id
  int <lower = 2> nlev_thresh; // Threshold for when we apply picking vs mult
  int wish_prior;  // 0 = Inverse Wishart, 1 = Wishart
  array[T] int<lower = 0, upper = 1> holdout; // 1=holdout, 0 =use
  int<lower = 0, upper = 1> check_bad; // 1 means to check bad
  array[T] int task_scale_group; // Always have at least 1 task_scale_group
}

transformed data{
  matrix[N, P] ind = rep_matrix(0, N, P); // ind_coded and ind_levels will expand to ind (below)
  matrix[P_cov, I] i_covt = i_cov'; // transpose of i_cov is what we use
  int df_tot = P + df;
  matrix[P,P] L_IW = cholesky_decompose(prior_cov * prior_cov_scale) * sqrt(df-1); // For Inverse Wishart
  matrix[P,P] L =    cholesky_decompose(prior_cov * prior_cov_scale) * 1/sqrt(df_tot); // For Wishart

  // Bartlett stuff
  int tri_n = 0;
  int P_IW = P;
  int P_W = 0;
  if (wish_prior == 1){
    P_IW = 0;
    P_W = P;
    tri_n = (P_W * (P_W - 1)) %/% 2; // # of lower tri elements 
  }
  array[P_W] real df_chi;
  array[tri_n,2] int tri_pos = get_pos_bart(P_W, P_W, tri_n); // position of bart_z elements

  int con_n_pos = vec_pos(con_sign);
  int con_n_neg = vec_neg(con_sign);
  int con_n = con_n_pos + con_n_neg;
  
  array[con_n]     int con_p;       // Array of which parameters are sign constrained
  array[con_n_pos] int con_p_pos;  // Positive constrained
  array[con_n_neg] int con_p_neg;  // Negative constrained
  matrix[con_n, I] con_delta;     // Con function scale for parameter and respondent
  int paircon_use = 0;
  vector[P] a_LB = rep_vector(-20, P); // LB for alpha (initial)
  vector[P] a_UB = rep_vector(20, P); // UB for alpha (initial)
  vector[P] prior_alpha_mod = prior_alpha; // Will modify alpha for constraints

  array[I] int resp_slice;  // Parallel threads slice across respondents 1:T  
  int count = 1; int count_pos =1; int count_neg = 1;
  array[I] real ll0_resp; // Compute LL of util = 0 for each respondent
  for (resp in 1:I) {
    real ll0_cum = 0;
    for (t in id_ranges[resp, 1]:id_ranges[resp, 2]){
        if (holdout[t] == 0) {
          ll0_cum += -log(end[t]-start[t]+1) * sum(dep[start[t]:end[t]]) * wts[t]; // log(1/nalt) = -log(nalt)
        }
    } // end t loop through tasks
    ll0_resp[resp] = ll0_cum;  // LL resp across tasks  
  }
  for (r in 1:I){resp_slice[r] = r;} // Slice by respondents
  for (i in 1:P){
    if (abs(con_sign[i]) > .00000001){ 
      con_p[count] = i; // array of parameters that are constrained
      con_delta[count,1:I] = rep_row_vector(con_sign[i] * con_factors[1], I); // con_delta controls hinge function
      count += 1;
      if (con_sign[i] > .00000001) {con_p_pos[count_pos] = i; count_pos +=1;}
      if (con_sign[i] < -.00000001) {con_p_neg[count_neg] = i; count_neg += 1;}
    }
    if (wish_prior == 1) df_chi[i] = df_tot - i + 1; // DF for prior on Bartlett diagonal 
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

  // Data transformations needed for speedup: multiplication vs picking
  // Indicators for which attributes and columns are used for picking (1 vs mult = 0)  
  array[n_atts] int atts_mult_bool = rep_array(1,n_atts); // How many attributes we do picking instead of multiplication
  array[P] int      cols_mult_bool = rep_array(1, P);  // 0/1 indicate which P we use for multiplication (vs Picking = 0)
  for (i in 1:n_atts){
    if (((code_blocks[i,2] - code_blocks[i,1]) > nlev_thresh) && code_blocks[i, 5] == 0){
      atts_mult_bool[i] = 0;  // att picked rather than multiply
      cols_mult_bool[code_blocks[i,1]:code_blocks[i,2]] = rep_array(0,code_blocks[i,2] - code_blocks[i,1] +1) ; // cols picked rather than multiply
    } 
  }
  array[sum(cols_mult_bool)] int cols_mult_pos; // Convert boolean cols to array of positions used for mult 
  int pos_index = 1;
  for (p in 1:P){
    if (cols_mult_bool[p] == 1){
      cols_mult_pos[pos_index] = p; // position index
      pos_index += 1; 
    } 
  }

  // Coding of attributes with levels
  // New array beta_full_pick stores which beta positions are picked in back coded beta
  int n_atts_pick = n_atts - sum(atts_mult_bool); // How many atts are picked
  array[N, n_atts_pick] int beta_full_pick;  // This will store the beta positions that are picked (vs mult) 
  matrix[sizes[3], P] code_master_block = rep_matrix(0, sizes[3], P);
  
  int lev_col = 1;  // Counter for which categorical attribute we are evaluating
  int att_pick = 1; // New counter for which picking attribute we are evaluating 
  for (i in 1:n_atts){
    int col_beg = code_blocks[i, 1];
    int col_end = code_blocks[i, 2];
    int row_beg = code_blocks[i, 3];
    int row_end = code_blocks[i, 4];
    int col_coded = code_blocks[i, 5]; // beginning of coded column
    
    code_master_block[row_beg:row_end,col_beg:col_end] = code_master[row_beg:row_end,col_beg:col_end]; 
    // code_master_block: each attribute defined by ONLY its own block of rows and cols, not interactions

    if (col_coded > 0){ //already coded, just copy
      ind[:,col_beg:col_end] = ind[:,col_beg:col_end] + ind_coded[:,col_coded:(col_coded + col_end - col_beg)];
    } else {
      int nrows = row_end - row_beg +1;
      matrix[nrows, P] code_att =  code_master[row_beg:row_end, ]; // code for this attribute
      for (n in 1:N){
        int lev = ind_levels[n,lev_col];
        if (lev >= 1) ind[n,] = ind[n,] + code_att[lev,];
      } // end 1:N loop coding        
      if (atts_mult_bool[i] == 0) { // Additional coding of beta_full_pick if picking
        for (n in 1:N){
          if (ind_levels[n,lev_col] == 0) {beta_full_pick[n,att_pick] = 1;
          } else {beta_full_pick[n,att_pick] = ind_levels[n,lev_col] + row_beg;}  // Absolute position + 1 in full recoded betas
          // Full back coded betas will have 1 additional level, with first coef = 0
        } // end 1:N loop coding
        att_pick += 1; // Next attribute for picking
      }
      lev_col = lev_col + 1; // Next categorcial level to code
    } // end if for coding att  
  } // end for (att coding)
   matrix[N, sum(cols_mult_bool)] ind_mult = ind[,cols_mult_pos]; // Remove columns used in picking
   print("Positive Constraints in Coded Cols: ", con_p_pos);
   print("Negative Constraints in Coded Cols: ", con_p_neg);
   print("Computing of each attribute (pick = 0, mult =1): ", atts_mult_bool);
   int n_scales = max(task_scale_group); 
} // end transformed data


parameters {
  vector[P] alpha; // upper mean
  matrix[P, P_cov] i_cov_load; // loadings for i_cov (resp covariates)
  matrix[P, I] z; // individual z scores N(0,1)

  // Inverse Wishart
  cholesky_factor_cov[P_IW] cov_chol_IW;

  // Wishart
  array[P_W] real<lower=0> bart_c ; // Bartlett diag^2
  array[tri_n] real bart_z; // Bartlett lower tri
  
  // Bad respondents
  array[check_bad] real<lower = 0, upper = 1> prob_bad;  // % of bad respondents
  array[n_scales-1] real<lower=0.001> scale_factor; // scale_factor for task_scale_group 1:n_scales
}


model {
  real LL_data;
  matrix[P,P] cov_chol;
  if (n_scales > 1) scale_factor ~ normal(1,1);
  
  if (wish_prior == 1) {
    if (tri_n == 0){ // No off-diagonal
        cov_chol = diag_post_multiply(L, to_vector(sqrt(bart_c))); // L * diag_matrix()
    } else {
        cov_chol = L * make_chol(sqrt(bart_c),bart_z,tri_pos); 
    }
    target+= chi_square_lupdf(bart_c|df_chi); // for (i in 1:P) bart_c[i] ~ chi_square(P + df - i + 1);
    if (tri_n > 0) bart_z ~ std_normal();
  } else {
    cov_chol = cov_chol_IW;
    target+= inv_wishart_cholesky_lupdf(cov_chol|df_tot, L_IW);
  }
  
  target+= normal_lupdf(alpha|prior_alpha_mod, a_sig);
  target += -sum(log1p_exp(-100 * (a_UB - alpha))); // Penalty if alpha > UB
  target += -sum(log1p_exp(-100 * (alpha - a_LB))); // Penalty if alpha < LB
  to_vector(z) ~ std_normal(); // respondent deviations from MVN(alpha, cov)

  if (P_cov > 0) to_vector(i_cov_load) ~ std_normal();
  if (check_bad == 1) target += normal_lupdf(prob_bad | .1, .3); 

  if (splitsize_n == 1){
    LL_data = reduce_sum(MNL_LL_resp_slice, resp_slice, splitsize[1],
                       ind_mult, dep, wts, holdout, start, end, npos_m1,
                       task_individual, id_ranges,
                       code_master_block, n_atts_pick, beta_full_pick, cols_mult_pos, // for pick atts
                       alpha, i_cov_load, i_covt, cov_chol, z,
                       paircon_use, con_n, con_p, con_delta, paircon_matrix,
                       con_p_pos, con_p_neg, con_factors,
                       check_bad, prob_bad, ll0_resp,scale_factor,task_scale_group);
  } else {
    profile("1"){
    LL_data = reduce_sum(MNL_LL_resp_slice, resp_slice, splitsize[1], // make sure matches - no way to do for loop
                       ind_mult, dep, wts, holdout, start, end, npos_m1,
                       task_individual, id_ranges,
                       code_master_block, n_atts_pick, beta_full_pick, cols_mult_pos, // for pick atts
                       alpha, i_cov_load, i_covt, cov_chol, z,
                       paircon_use, con_n, con_p, con_delta, paircon_matrix,
                       con_p_pos, con_p_neg, con_factors,
                       check_bad, prob_bad, ll0_resp,scale_factor,task_scale_group);
    }
    profile("2"){
    LL_data = reduce_sum(MNL_LL_resp_slice, resp_slice, splitsize[2], // make sure matches - no way to do for loop
                       ind_mult, dep, wts, holdout, start, end, npos_m1,
                       task_individual, id_ranges,
                       code_master_block, n_atts_pick, beta_full_pick, cols_mult_pos, // for pick atts
                       alpha, i_cov_load, i_covt, cov_chol, z,
                       paircon_use, con_n, con_p, con_delta, paircon_matrix,
                       con_p_pos, con_p_neg, con_factors,
                       check_bad, prob_bad, ll0_resp,scale_factor,task_scale_group);
    }
    if (splitsize_n >= 3) {
      profile("3"){
      LL_data = reduce_sum(MNL_LL_resp_slice, resp_slice, splitsize[3], // make sure matches - no way to do for loop
                       ind_mult, dep, wts, holdout, start, end, npos_m1,
                       task_individual, id_ranges,
                       code_master_block, n_atts_pick, beta_full_pick, cols_mult_pos, // for pick atts
                       alpha, i_cov_load, i_covt, cov_chol, z,
                       paircon_use, con_n, con_p, con_delta, paircon_matrix,
                       con_p_pos, con_p_neg, con_factors,
                       check_bad, prob_bad, ll0_resp,scale_factor,task_scale_group);
    }}
    if (splitsize_n >= 4) {
      profile("4"){
      LL_data = reduce_sum(MNL_LL_resp_slice, resp_slice, splitsize[4], // make sure matches - no way to do for loop
                       ind_mult, dep, wts, holdout, start, end, npos_m1,
                       task_individual, id_ranges,
                       code_master_block, n_atts_pick, beta_full_pick, cols_mult_pos, // for pick atts
                       alpha, i_cov_load, i_covt, cov_chol, z,
                       paircon_use, con_n, con_p, con_delta, paircon_matrix,
                       con_p_pos, con_p_neg, con_factors,
                       check_bad, prob_bad, ll0_resp,scale_factor,task_scale_group);
    }}
    if (splitsize_n >= 5) {
      profile("5"){
      LL_data = reduce_sum(MNL_LL_resp_slice, resp_slice, splitsize[5], // make sure matches - no way to do for loop
                       ind_mult, dep, wts, holdout, start, end, npos_m1,
                       task_individual, id_ranges,
                       code_master_block, n_atts_pick, beta_full_pick, cols_mult_pos, // for pick atts
                       alpha, i_cov_load, i_covt, cov_chol, z,
                       paircon_use, con_n, con_p, con_delta, paircon_matrix,
                       con_p_pos, con_p_neg, con_factors,
                       check_bad, prob_bad, ll0_resp,scale_factor,task_scale_group);
    }}    
  } // end splitsize test
  target += LL_data;
} // End Model

generated quantities{
  matrix[P, I] beta_ind;
  {
  matrix[P,P] cov_chol;
  if (wish_prior == 1) {
    if (tri_n == 0){ // No off-diagonal
        cov_chol = diag_post_multiply(L, to_vector(sqrt(bart_c))); // L * diag_matrix()
    } else {
        cov_chol = L * make_chol(sqrt(bart_c),bart_z,tri_pos); 
    }
  } else {cov_chol = cov_chol_IW;}
  beta_ind = rep_matrix(alpha, I) + (i_cov_load * i_covt) + cov_chol * z; // Unconstrained
    if (con_n > 0){
      beta_ind[con_p,1:I] = con_delta .* log1p_exp(beta_ind[con_p,1:I] ./ con_delta);
    } 
  }
}

