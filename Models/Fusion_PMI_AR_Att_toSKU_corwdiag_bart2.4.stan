// Kevin Lattery Jan 2022
// Store Level sales data (Over time) + Conjoint Data
// Uses sku attribute decomposition d_att

functions{
  matrix logistic_hinge(matrix x, matrix delta) {
    return delta .* log1p_exp(x ./ delta);
    // https://statmodeling.stat.columbia.edu/2017/05/19/continuous-hinge-function-bayesian-modeling/
  } // Not used here
  

  // Function used for initial AR delivered June 2021 and Dec 2021
  real AR_MNLwSimPop(int[] array_slice,
                   int a_beg, int a_end, // Stan determines values of these for specific core
                   matrix sim_pop,
                   vector lag_share, vector u_sku_from_att, int[] d_sku,
                   vector b_price, vector lag_price, vector new_price,
                   vector b_dist, vector lag_dist, vector new_dist,
                   vector b_trend, vector time_diff,
                   vector dep, vector dep_att,
                   int[] start, int[] end, int[] end_lag,
                   real ar_scale
  ) {
    int npop = cols(sim_pop); // each column (not row) is a respondent
    real loglike = 0;
    for (t in a_beg:a_end){
     int nalt_task = (end[t] - start[t] + 1); // # of alternatives in new/forecast task
     int sku_t[nalt_task] = d_sku[start[t]:end[t]]; // skus in new/forecast task
     
     int nalt_task_lag = (end_lag[t] - start[t] + 1); // # of alternatives in lag task
     int sku_t_lag[nalt_task_lag] = d_sku[start[t]:end_lag[t]]; // lag skus in task
     vector[nalt_task_lag] lag_share_log = log(lag_share[start[t]:end_lag[t]]); // lag log shares
     
      vector[nalt_task_lag] u_global_lag = 
                  (b_price[sku_t_lag] .* lag_price[start[t]:end_lag[t]]) +
                  (b_dist[sku_t_lag]  .* lag_dist[start[t]:end_lag[t]]);
                  
      vector[nalt_task] u_global_new = 
                  (b_price[sku_t] .* new_price[start[t]:end[t]]) +
                  (b_dist[sku_t]  .* new_dist[start[t]:end[t]]);
                  
      vector[nalt_task] u_trend = b_trend[sku_t] .* time_diff[start[t]:end[t]];
                  
      vector[nalt_task] pred_share = softmax(u_sku_from_att[sku_t] + u_global_new); // predicted forecast just from attributes
      vector[nalt_task] mu_new; // fill in below

     // Compute shares for task t for each member of population
     vector[nalt_task] task_prob_sum_new = rep_vector(0,nalt_task);
     vector[nalt_task_lag] task_prob_sum_lag = rep_vector(0,nalt_task_lag);
     vector[nalt_task] U_final;
     vector[nalt_task] pred_new;
     
     // Define global mean for new/forecast period
     mu_new = log(pred_share); // set overlap skus below
     mu_new[1:nalt_task_lag] = lag_share_log + log(sum(pred_share[1:nalt_task_lag])) +
                               (u_global_new + u_trend)[1:nalt_task_lag] - u_global_lag;
     // do simulations with mu = lag_share_log, mu_new
     for (i in 1:npop){
       task_prob_sum_new += softmax((col(sim_pop,i)[sku_t]     + mu_new)/ar_scale);
       task_prob_sum_lag += softmax((col(sim_pop,i)[sku_t_lag] + lag_share_log)/ar_scale);
     }
     // Compute final Utility and LogLikelihood
     pred_new = task_prob_sum_new/npop;
     U_final = log(pred_new);
     U_final[1:nalt_task_lag] = lag_share_log + log(sum(pred_new[1:nalt_task_lag])) +
                                log(pred_new[1:nalt_task_lag]) - log(task_prob_sum_lag/npop);
     loglike += dot_product(log_softmax(U_final), dep[start[t]:end[t]]); // VAR model LL of task
     loglike += dot_product(log(pred_share),dep_att[start[t]:end[t]]); // Simple attribute LL (sep wts)
    }
    return loglike;
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
  // 1 = Real World, 2 = Conjoint, Both are units
  int N; // # Rows in stacked data file
  int P; // # Coded attribues
  int P_sku;  // # SKUs

  // Main data
  // Real World Data
  int<lower = 1, upper = P_sku> d_sku[N]; // SKUs (array)
  vector[N] lag_share; 
  vector[N] lag_price;
  vector[N] lag_dist;
  vector[N] new_price;
  vector[N] new_dist;
  vector[N] time_diff;
  vector<lower = 0, upper = 1>[N] dep; // Dep variable
  vector[N] wts; // weights of each row
  vector[N] wts_att; // weights for LL fit of simple attribute model
  matrix[P_sku, P] sku_to_att;   // converts skus to attributes
  cov_matrix[P] prior_cor;  // Prior covariance/correlation matrix from conjoint
  
  // Ragged array matching, For each task, Specify:
  int T; // # Tasks
  int<lower = 1, upper = N> start[T]; // which row task begins
  int<lower = 1, upper = N> end[T];   // which row task ends
  int<lower = 1, upper = N> end_lag[T];   // which row the lag task ends

  // Raw simulated data ~ N(0,1)
  int I; // # How many simulated buyers we want for estimation
  matrix[P, I] sim_pop_z; // Random members of pop ~ N(0,1), converts to sim pop
  
  // Sigma priors 
  real<lower = 0> price_sigma;  // b_price ~ N(p_mu, price_sigma) where p_mu is parameter
  real p_mu_mu; // p_mu ~ N(p_mu_mu,p_mu_sigma)
  real<lower = 0> p_mu_sigma; 
  real<lower = 0> trend_sigma;  // b_trend ~ N(0, trend_sigma)
  real<lower = 0> dist_sigma;  // b_dist ~ N(0, dist_sigma)
  real<lower = 0> cov_diag_sigma; // cov_diag ~ N(,cov_diag_sigma)
  
  real<lower = .25, upper = 1> ar_scale; // Scaling of Smoothed Accept-Reject
  int<lower = 0> df;
  
  // Other Data
  matrix[P,P] cov_block; // Specifies blocks of covariance items
  int<lower = 0> splitsize; // grainsize for parallel processing
}

transformed data{
  vector[N] dep_wt = dep .* wts;
  vector[N] dep_wt_att = dep .* wts_att;
  matrix[P,P] L = cholesky_decompose(prior_cor/(P + df)); // chol(cor)/sqrt(DF) 
  int tri_n = tri_sum(cov_block); // # of lower tri elements 
    int tri_pos[tri_n,2] = get_pos(cov_block, tri_n); // position of elements
  
  int array_slice[T]; // Parallel threads slice (stores)  
  real df_chi[P];
  for (i in 1:P){
    df_chi[i] = P + df - i + 1;
  }  
  for (i in 1:T){
    array_slice[i] = i;
  }
}

parameters {
  // Real World model
  vector<upper = -.5>[P_sku] b_price; // price utilities each item/sku (agg level only)
  vector[P_sku] b_trend; // Each sku will have a trend slope
  vector<lower = .5, upper = 2>[P_sku] b_dist; // distribution coefficient
  real<upper = -.5> p_mu; // mean of price beta
  
  // Cov matrix
  // vector<lower=0>[P] cov_diag; // diag of cov
    // Cov matrix
  vector<lower=0>[P] cov_diag; // diag of cov
  real bart_z [tri_n]; // Bartlett lower tri
  real<lower=0> bart_c [P]; // Bartlett diag^2

  vector<lower = .5, upper = 5>[P] sd_diag; // separate diagonal
  
  vector[P] b_attributes; // coefficients to predict share without AR
}

transformed parameters {
  vector[P_sku] u_sku_from_att = (sku_to_att * b_attributes); // Utility of sku based on b_attributes
  matrix[P,P] cov_chol = diag_pre_multiply(sd_diag, L) * make_chol(sqrt(bart_c),bart_z,tri_pos);
}

model {
  // create simulated population of skus
  matrix[P_sku, I] sim_pop_raw = sku_to_att * (cov_chol * sim_pop_z); // sim_pop from atts converted to skus
  vector[P_sku] mu_pop = inv(I) * (sim_pop_raw * rep_vector(1, I)); // sku means
  matrix[P_sku, I] sim_pop = sim_pop_raw - rep_matrix(mu_pop, I) ; // 0 centered sku ppopulation
  // matrix[P,P] prior_scale_inv = quad_form_diag(prior_cor_inv, inv(sd_diag)) * (P + df);
  // inv(scale = cov/k) = inv((D * Corr *D)/k)  = (1/D * Corr_inv * 1/D) * k, where k = P + df
  
  // priors on the parameters
  b_attributes ~ normal(0,10);
  b_price ~ normal(p_mu, price_sigma); // Truncated since declared bounded
  p_mu ~ normal(p_mu_mu,p_mu_sigma);
  b_trend ~ normal(0, trend_sigma);
  b_dist ~ normal(1, dist_sigma); 
  sd_diag ~ normal(1, 10);
  target+= chi_square_lpdf(bart_c|df_chi); // for (i in 1:P) bart_c[i] ~ chi_square(P + df - i + 1);
  to_vector(bart_z) ~ std_normal();

  //target += log_chol_wish(cov_chol, prior_scale_inv, df); # cov_chol ~ wish(cor, df)
  // target += log_chol_wish(cov_chol, prior_scale_inv, df) + (P+df) * log(determinant(prior_scale_inv));
  // since prior varies, subtract log[det(S) ^ (n/2)] = (-n/2) * log[det(Sinv)]

  // loglikelihood
  target += reduce_sum(AR_MNLwSimPop, array_slice, splitsize,
                   sim_pop,
                   lag_share, u_sku_from_att, d_sku,
                   b_price, lag_price, new_price,
                   b_dist, lag_dist, new_dist,
                   b_trend, time_diff,
                   dep_wt, dep_wt_att,
                   start, end, end_lag,
                   ar_scale
                   );
} // End Model


