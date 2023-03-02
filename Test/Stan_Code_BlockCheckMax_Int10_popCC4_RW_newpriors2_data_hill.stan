// Kevin Lattery Jan 2022
// Store Level sales data (Over time) + Conjoint Data
// Uses sku attribute decomposition d_att

functions{
  matrix logistic_hinge(matrix x, matrix delta) {
    return delta .* log1p_exp(x ./ delta);
    // https://statmodeling.stat.columbia.edu/2017/05/19/continuous-hinge-function-bayesian-modeling/
  } // Not used here
  
  vector row_max(matrix x){
    vector[rows(x)] result;
    for (i in 1:rows(x)){
      result[i] = max(x[i,]);      
    }
    return result;
  }
  
  vector col_max(matrix x){
    vector[cols(x)] result;
    for (i in 1:cols(x)){
      result[i] = max(x[,i]);      
    }
    return result;
  }
  
  real log_chol_wish(matrix chol_x, matrix prior_scale_inv, int df_overp) {
    // chol_x is    p x p cholesky of matrix x to test 
    // prior_scale_inv is p x p inverse of prior scale matrix
    //   I use prior_scale_inv = inv(prior_cov/total_df)
    // df_overp is 1+ degrees of freedom over p (total df = p + df_overp)

    // returns log(wishart(x | prior_inv, p + df_overp))
    //    excluding pdf normalizing constants from df_overp, prior_scale_inv
     real log_det_x = sum(log(square(diagonal(chol_x)))); // log(determinant(x))
     real log_e_num = -.5 * trace(prior_scale_inv * tcrossprod(chol_x)); // log(e^numerator) = numerator
     return ((df_overp - 1.0) * .5  * log_det_x + log_e_num);
   }
   
   vector log_w0(vector index_binary, vector x){
     return(index_binary .* log(x + 1.0000000001 - index_binary));
     // index = 1 returns log(x), index = 0 returns 0
   }
   vector exp_lam_w0(vector index_binary, vector x, real lam){
     return(index_binary .* exp(lam * log(x + 1.0000000001 - index_binary)));
     // index = 1 returns log(x), index = 0 returns 0
   }
   
   vector maxexp_u_whatif(vector sim_pop_one,
                    vector mu_self, vector mu_whatif,
                    matrix morph_rowtocol, vector morph_receive,
                    vector mult_self_b, vector mult_from_b, int P_sku){
          vector[P_sku] exp_util =  mult_self_b .* exp(mu_self + sim_pop_one);
          vector[P_sku] exp_util_whatif = exp(mu_whatif + sim_pop_one);
          for (i in 1:P_sku){
             if (morph_receive[i] > .99){
                 exp_util[i] = fmax(exp_util[i],
                             mult_from_b[i] * max(morph_rowtocol[,i] .* exp_util_whatif) +
                         (1-mult_from_b[i]) * exp_util[i]);
             }
          }
          return(exp_util);
    }
    
  // Created May 2022 to allow morphing
  real AR_MNLwSimPop_max(int[] array_slice,
                   int a_beg, int a_end, // Stan determines values of these for specific core
                   matrix sim_pop_int, row_vector sim_pop_price, matrix[] morph_rowtocol,
                   vector[,] share, 
                   vector[,] log_price,
                   vector b_dist, vector[,] log_dist,
                   vector b_trend,
                   vector int_sku_over, vector int_sku_npl,
                   real b_promo_exist, real b_promo_npl, vector[,] promo,
                   real b_aware, vector[,] log_aware,
                   real b_morph_npl_self, real b_morph_npl_hyp, real b_morph_exist_hyp,
                   vector b_channel_exist, vector b_channel_npl, matrix[] channel_bin,
               //    vector b_tier_sku, vector[,] tier_price_sku,
                   int per_base, int[] per_lag, int[] per_new, int[] region,
                   vector[,] skus_bin, vector lag_use,
                   vector[,] wts, vector[,]wts_att,
                   real hill_g, real hill_a // hill function
  ) {
    int npop = cols(sim_pop_int); // each column (not row) is a respondent
    // int P = rows(sim_pop_int);
    int P_sku = rows(sim_pop_int);
    real loglike = 0;
    for (t in a_beg:a_end){
        int p_lag = per_lag[t]; // lag period for task
        int p_new = per_new[t]; // new/forecast period for task
        int region_t = region[t];
        vector[P_sku] skus_new = skus_bin[region_t, p_new]; // skus in new/forecast period: absolute
        //vector[P_sku] skus_lag = skus_bin[region_t, p_lag]; // skus in lag period: absolute
        vector[P_sku] skus_lag = lag_use .* floor(
          (skus_bin[region_t, p_lag] + skus_bin[region_t, p_lag-1] + skus_bin[region_t, p_lag-2] + .1)/3
          );
        
        vector[P_sku] skus_over = skus_new .* skus_lag; // skus in both periods: relative
        vector[P_sku] skus_npl = skus_new - skus_over; // skus that are new product launches vs lag: relative
        vector[P_sku] mu_int = (skus_over .* int_sku_over) + (skus_npl .* int_sku_npl); // added
        vector[P_sku] b_promo = (skus_over * b_promo_exist) + (skus_npl * b_promo_npl); // added
        vector[P_sku] promo_hill = pow(promo[region_t, p_new], hill_a) ./
                                  (pow(promo[region_t, p_new], hill_a) + pow(hill_g, hill_a));           
        vector[P_sku] u_global_lag = mu_int + b_trend * (p_lag - per_base) + 
                  log_dist[region_t, p_lag] .* b_dist  +
               // tier_price_sku[region_t, p_lag] .* b_tier_sku +
                  promo_hill                 .* b_promo  +
                  log_aware[region_t, p_lag] * b_aware   +
                  (channel_bin[p_lag]        * b_channel_exist) .* skus_over +
                  (channel_bin[p_lag]        * b_channel_npl) .* skus_npl;
                  
         vector[P_sku] u_global_fore = mu_int + b_trend * (p_new - per_base) + 
                  log_dist[region_t, p_new] .* b_dist  +
               // tier_price_sku[region_t, p_new] .* b_tier_sku +                
                  promo_hill                 .* b_promo  +
                  log_aware[region_t, p_new] * b_aware   +
                  (channel_bin[p_new]        * b_channel_exist) .* skus_over +
                  (channel_bin[p_new]        * b_channel_npl) .* skus_npl;
          
        // vector[P_sku] u_trend = b_trend * (p_new - p_lag);
        // vector[P_sku] cal_lag_all = log_w0(skus_lag, share[region_t, p_lag]);
         vector[P_sku] share_lag_over = (share[region_t, p_lag] .* skus_over)/sum(share[region_t, p_lag] .* skus_over);
              ///////////////////

    // Simulations
     ////////////////////
   // Declare some storage
     vector[P_sku] task_prob_sum_new = rep_vector(0,P_sku);
     vector[P_sku] task_prob_sum_lag = rep_vector(0,P_sku);
     vector[P_sku] sim_share_fore;
     vector[P_sku] sim_share_lag;
     vector[P_sku] U_final;
     vector[P_sku] pred_new;
     vector[P_sku] exp_util; // placeholder for e^U calcs

   // Forecast
     matrix[P_sku,P_sku] morph_rowtocol_clean = diag_post_multiply(morph_rowtocol[p_new], skus_new);//  Remove any skus that no longer exist
     // util of sku @ price, dist, etc of sku morphing into
     // cannot use mu_int because products that no longer exist will have mu_int = att npl 
     // all what_if skus are not npl 
     vector[P_sku] u_global_whatif =  int_sku_over +
                      b_trend * (p_new - per_base) +
                    //  promo[region_t, p_new]     * b_promo  +
                    //  log_aware[region_t, p_new] * b_aware   +
                      channel_bin[p_new]         * b_channel_exist +
                   //   tier_price_sku[region_t, p_new] .* b_tier_sku +                       
                      (morph_rowtocol_clean * log_dist[region_t, p_new]) .* b_dist;
     vector[P_sku] mu_whatif = u_global_whatif .* row_max(morph_rowtocol_clean);

     vector[P_sku] morph_receive = col_max(morph_rowtocol_clean); // 0/1 skus receive                                      
     vector[P_sku] morph_npl = morph_receive .* skus_npl; // 0/1 NPLs that are morphed into
     vector[P_sku] morph_exist = morph_receive .* (1-skus_npl); // 0/1 existing skus morphed into
     vector[P_sku] mult_self_b = (skus_new - morph_npl) + morph_npl * b_morph_npl_self;
     vector[P_sku] mult_from_b = morph_npl * b_morph_npl_hyp + morph_exist * b_morph_exist_hyp;

     vector[P_sku] log_price_t = log_price[region_t, p_new];
     vector[P_sku] log_price_t_whatif = morph_rowtocol_clean * log_price_t;
     for (i in 1:npop){
        exp_util = skus_new .* maxexp_u_whatif(
                     sim_pop_int[,i],
                     u_global_fore + (sim_pop_price[i] * log_price_t), // change from v4
                     mu_whatif + (sim_pop_price[i] * log_price_t_whatif), // change from v4
                     morph_rowtocol_clean, morph_receive,
                     mult_self_b, mult_from_b, P_sku);
        task_prob_sum_new += exp_util/sum(exp_util); // key: sum of shares over pop            
     } 

     // Repeat for Lag  
    // Remove any skus not in overlap & that do not contribute to *FORECAST* (morph_self_clean forecast)
     morph_rowtocol_clean = diag_post_multiply(morph_rowtocol[p_lag], skus_over .* morph_receive); // yes this refers to forecast
     morph_receive = col_max(morph_rowtocol_clean); // Now redefine for lag
     // util of sku @ price, dist, etc of sku morphing into
     u_global_whatif =  int_sku_over +
               b_trend * (p_lag - per_base) +
              // promo[region_t, p_lag]     * b_promo  +
              //log_aware[region_t, p_lag] * b_aware   +
               channel_bin[p_lag]         * b_channel_exist +
            //   tier_price_sku[region_t, p_lag] .* b_tier_sku + 
               (morph_rowtocol_clean * log_dist[region_t, p_lag]) .* b_dist;
     mu_whatif = u_global_whatif .* row_max(morph_rowtocol_clean);

     morph_receive = col_max(morph_rowtocol_clean); // 0/1 skus receive                                      
     morph_npl = morph_receive .* skus_npl; // 0/1 NPLs that are morphed into
     morph_exist = morph_receive .* (1-skus_npl); // 0/1 existing skus morphed into
     mult_self_b = (skus_new - morph_npl) + morph_npl * b_morph_npl_self;
     mult_from_b = morph_npl * b_morph_npl_hyp + morph_exist * b_morph_exist_hyp;

     log_price_t = log_price[region_t, p_lag];
     log_price_t_whatif = morph_rowtocol_clean * log_price_t;
     for (i in 1:npop){
        exp_util = skus_over .* maxexp_u_whatif(
                     sim_pop_int[,i],
                     u_global_lag + (sim_pop_price[i] * log_price_t),
                     mu_whatif + (sim_pop_price[i] * log_price_t_whatif),
                     morph_rowtocol_clean, morph_receive,
                     mult_self_b, mult_from_b, P_sku);
         task_prob_sum_lag += exp_util/sum(exp_util); // key: sum of shares over pop 
     }
     
     // Compute final Utility and LogLikelihood
     sim_share_fore = task_prob_sum_new/npop;
     sim_share_lag = task_prob_sum_lag/npop;
     U_final = skus_over .* share_lag_over .* sim_share_fore ./ (sim_share_lag + 1 - skus_over) +
               skus_npl .* sim_share_fore; 
     pred_new = U_final/sum(U_final);
    // print(t, ";", 
    // dot_product(log_w0(skus_new, pred_new), wts[region_t, p_new] .* share[region_t, p_new]),
    // ";",
    // dot_product(log_w0(skus_new, sim_share_fore), wts_att[region_t, p_new] .* share[region_t, p_new]));
     
     loglike += dot_product(log_w0(skus_new, pred_new), wts[region_t, p_new] .* share[region_t, p_new]); // VAR model LL of task
     loglike += dot_product(log_w0(skus_new, sim_share_fore), wts_att[region_t, p_new] .* share[region_t, p_new]); // Non-Var

    } // end t
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
  int task_beg;
  int task_end;
  int P; // # Coded attributes
  int P_int; // Coded attributes related to sku intercept (not price)
  int P_sku;  // # SKUs
  int N_periods; // # Periods
  int N_regions;
  int<lower = 1> per_base;  // base period for intercept
  int<lower = 0> n_channels;
  int<lower = 0> nbrand; // number of brand levels

  // Main data
  // Real World Data
  vector[P_sku] price[N_regions, N_periods];
  vector[P_sku] dist[N_regions, N_periods];
  vector[P_sku] aware[N_regions, N_periods];
  vector[P_sku] promo[N_regions, N_periods];
  vector[P_sku] share[N_regions, N_periods];
  vector[P_sku] wts[N_regions, N_periods];
  vector[P_sku] wts_att[N_regions, N_periods];
  vector<lower = 0, upper = 1> [P_sku] skus_bin[N_regions, N_periods];  // 1 = include in forecast, 0 = exclude
  matrix<lower = 0, upper = 1> [P_sku, n_channels] channel_bin[N_periods];

  matrix[P_sku, P_int] sku_to_att;   // converts skus to attributes, Pth paramter is price
  cov_matrix[P] prior_cor;  // Prior covariance/correlation matrix from conjoint

  // Raw simulated data ~ N(0,1)
  int I; // # How many simulated buyers we want for estimation
  matrix[P, I] sim_pop_z; // Random members of pop ~ N(0,1), converts to sim pop
  
  matrix[P_sku,P_sku] morph_rowtocol[N_periods]; // Array of matrices for each period

  int T; // Number of tasks
  int<lower = 1, upper = N_periods> per_lag[T]; // lag period number for each task
  int<lower = 1, upper = N_periods> per_new[T]; // new
  int<lower = 1, upper = N_regions> region[T]; // 

  // Sigma priors 
  real<lower = 0> price_sigma;  // b_price ~ N(p_mu, price_sigma) where p_mu is parameter
  real p_mu_mu; // p_mu ~ N(p_mu_mu,p_mu_sigma)
  real<lower = 0> p_mu_sigma; 
  real<lower = 0> trend_sigma;  // b_trend ~ N(0, trend_sigma)
  real<lower = 0> dist_sigma[P_sku];  // b_dist ~ N(0, dist_sigma)
  real<lower = 0> cov_diag_sigma; // cov_diag ~ N(,cov_diag_sigma)
  real<lower = 0> att_sigma[P_int]; // Prior on attribute shrinkage, exclude Price
  real<lower = 0> b_promo_prior_mu[2]; // b_promo_exist/np; ~ N(mu, 1): .1, .2 Bristol, .5 RW
  
  real<lower = .25, upper = 1> ar_scale; // Scaling of Smoothed Accept-Reject
  int<lower = 0> df;
  real<lower = 1> lam_inv; // nested logit 1/lam
  
  // Other Data
  matrix[P,P] cov_block; // Specifies blocks of covariance items
  int<lower = 0> splitsize; // grainsize for parallel processing
  vector<lower = 0, upper = 1>[P_sku] lag_use; // 0 = do not use sku lag, 1 = use as normal
  
//  int<lower = 1> n_tiers;  // Price Tier x Format or other groups
//  int<lower = 1, upper = n_tiers> sku_tier[P_sku]; // Classification of sku into tier
//  vector[n_tiers] tier_price[N_regions, N_periods]; // Mean price of each tier
//  vector[P_sku] price_rel[N_regions, N_periods]; // this will be relative price in tier now
}

transformed data{
  vector[P_sku] b_trend = rep_vector(0, P_sku);
  vector[2] sli_minus_reg = [2,1]'; // sli [1,0] minus reg [-1,-1] 
  real size100_minus_ks = 2; // [1] - [-1]
  vector[P_sku] log_price[N_regions, N_periods];
  vector[P_sku] log_dist[N_regions, N_periods];
  vector[P_sku] log_share[N_regions, N_periods];
  vector[P_sku] log_aware[N_regions, N_periods];
  vector[P_sku] dep_wt[N_regions, N_periods];
  vector[P_sku] dep_wt_att[N_regions, N_periods];
  vector[P_sku] tier_price_sku[N_regions, N_periods];
  
  matrix[P,P] L = cholesky_decompose(prior_cor/(P + df)); // chol(cor)/sqrt(DF) 
  int tri_n = tri_sum(cov_block); // # of lower tri elements 
    int tri_pos[tri_n,2] = get_pos(cov_block, tri_n); // position of elements
  
  int array_slice[task_end - task_beg + 1]; // Parallel threads slice (tasks)  
  real df_chi[P];
  for (i in 1:N_regions){
    for (j in 1:N_periods){
      log_price[i,j] = log_w0(skus_bin[i,j], price[i,j] * .01); // log, but 0 when skus_bin = 0
      log_dist[i,j] = log_w0(skus_bin[i,j], dist[i,j] * 10);
      log_share[i,j] = log_w0(skus_bin[i,j], share[i,j]);
      log_aware[i,j] =  log(aware[i, j] + .05) - log(1.05);
      dep_wt[i,j] = wts[i,j] .* share[i,j];
      dep_wt_att[i,j] = wts_att[i,j] .* share[i,j];
    //  tier_price_sku[i,j] = tier_price[i,j][sku_tier]; // mean prices @tier level to sku level
    }
  }
  for (i in 1:P){
    df_chi[i] = P + df - i + 1;
  }  
  for (i in 1:(task_end - task_beg + 1)){ // SHOULD BE T
    array_slice[i] = task_beg + i - 1;
  }
}

parameters {
  // Real World model
//  vector[P_sku] b_trend; // Each sku will have a trend slope
  vector<lower = .5, upper = 1.5>[P_sku] b_dist; // distribution coefficient
  real<upper = -.5> p_mu; // mean of price beta
  
  real<lower = 0> b_promo_exist; // beta for sku_promo ln(dist * e^b_promo)
  real<lower = 0> b_promo_npl; // beta for sku_promo ln(dist * e^b_promo)

  real<lower = .01, upper = 1> b_aware; // convert aware% to a log scale: ln(1/b_aware+aware)/ln(1/b_aware + 1)

  real<lower = .75, upper = 1>  b_morph_npl_self; // % of full npl
  real<lower =.25, upper = 1> b_morph_npl_hyp; // % from product hypothetical A, new B
  real<lower =0, upper = .5>  b_morph_exist_hyp; // % from hypthetical A, existing B 
  
  // Cov matrix
  // vector<lower=0>[P] cov_diag; // diag of cov
    // Cov matrix
  // vector<lower=.5, upper = 2>[P] cov_diag; // diag of cov
  real bart_z [tri_n]; // Bartlett lower tri
  real<lower=0> bart_c [P]; // Bartlett diag^2

  vector<lower = .5, upper = 4>[P] sd_diag; // separate diagonal was .5, 5
  
  vector[P_int] b_attributes_over; // coefficients to get mean of utilities for overlap
  vector[8] b_attributes_nplspec;  // coefficients specific to npl 
  real<lower = .1, upper = 1> k_brand; // scale factor to adjust brand <= 1
  real<lower = .1, upper = 1> k_fruit_cps;
  real<lower = .1, upper = 1> k_flavor;
  real<lower = .1, upper = 1> k_color;
  vector[n_channels] b_channel_exist; // 1 beta per channel
  
  real<lower = 1, upper = 3> hill_a; // hill function exponent
  real<lower = .01, upper = 1> hill_g; // hill_function gamma
  
 // vector<upper = 0> [n_tiers] b_tier; // price elasticity of tier (negative)
}

transformed parameters {
  vector[n_channels] b_channel_npl =  b_channel_exist; // 1 beta per channel
  matrix[P,P] cov_chol = diag_pre_multiply(sd_diag, L) * make_chol(sqrt(bart_c),bart_z,tri_pos);

  vector[P_sku] int_sku_over = (sku_to_att * b_attributes_over); // Utility of sku based on b_attributes
  
  vector[P_sku] int_sku_npl;
  vector[P_int] b_attributes_npl;
  // vector[P_sku] b_tier_sku = b_tier[sku_tier]; // tier parameter for each sku
  
  // This section customized for each data set
  b_attributes_npl[1:nbrand] = k_brand * b_attributes_over[1:nbrand]; // scale brand attributes
  b_attributes_npl[nbrand+1] = b_attributes_nplspec[1]; // Length
  b_attributes_npl[(nbrand+2):(nbrand+4)] = b_attributes_nplspec[2:4]; // Tier1-3
  b_attributes_npl[(nbrand+5):(nbrand+6)] = b_attributes_nplspec[5:6]; // Tier1_Format
  b_attributes_npl[(nbrand+7):(nbrand+8)] = b_attributes_nplspec[5:6]; // Tier2_Format
  b_attributes_npl[(nbrand+9):(nbrand+10)] = b_attributes_nplspec[5:6]; // Tier3_Format
  b_attributes_npl[(nbrand+11):(nbrand+12)] = b_attributes_nplspec[7:8]; // Tier4_Format
  b_attributes_npl[(nbrand+13):(nbrand+14)] = k_fruit_cps * b_attributes_over[(nbrand+13):(nbrand+14)]; // fruit_cps
  b_attributes_npl[(nbrand+15):(nbrand+20)] = k_flavor * b_attributes_over[(nbrand+15):(nbrand+20)];  // Flavor
  b_attributes_npl[(nbrand+21):(nbrand+23)] = k_color * b_attributes_over[(nbrand+21):(nbrand+23)];  // Color
  
  int_sku_npl = (sku_to_att * b_attributes_npl);
}

model {
  // create simulated population of skus
  matrix[P, I] sim_pop_att = cov_chol * sim_pop_z; // sim_pop on attributes
  row_vector[I] sim_pop_price = -.2 * log1p_exp(inv(-.2) * (p_mu + sim_pop_att[P,])); // Pop price utilities constrained
  
  matrix[P_sku, I] sim_pop_int = sku_to_att * sim_pop_att[1:P_int,]; // Pop sku utilities from attributes

  // priors on the parameters
  b_attributes_over ~ normal(0,att_sigma);
  b_attributes_nplspec ~ normal(0,10); 
  p_mu ~ normal(p_mu_mu, p_mu_sigma);
//  b_trend ~ normal(0, trend_sigma);
  b_dist ~ normal(1, dist_sigma); 
  sd_diag ~ normal(1, 2); // was 1,5
  b_promo_exist ~ normal(b_promo_prior_mu[1],1);    
  b_promo_npl ~ normal(b_promo_prior_mu[2],1);   
  b_aware ~ normal(.5,1);
  to_vector(b_channel_exist) ~ normal(0,5);
  to_vector(b_channel_npl) ~ normal(0,5);
  target+= chi_square_lpdf(bart_c|df_chi); // for (i in 1:P) bart_c[i] ~ chi_square(P + df - i + 1);
  to_vector(bart_z) ~ std_normal();
  b_morph_npl_self ~ beta(1.2, 1); // scale factor discount from npl 
  b_morph_npl_hyp ~ beta(1.1, 1.1); // % from product A
  b_morph_exist_hyp ~ beta(1.2, 1.3); // mode = .4 
  k_brand ~ normal(1,5); // scale factor to adjust brand
  k_fruit_cps ~ normal(1,5); // scale factor to adjust brand
  k_color ~ normal(1,5); // scale factor to adjust brand
  k_flavor ~ normal(1,5);
  hill_a ~ normal(2,1); // hill function priors
  hill_g ~ beta(1.05,1.05); // Quite vague
  
//  b_tier ~ normal(-1, 5); // price elasticity of tiers
  
  target += -(log1p_exp(-100 * dot_product(sli_minus_reg, b_attributes_nplspec[5:6])));
  target += -(log1p_exp(-100 * dot_product(sli_minus_reg, b_attributes_nplspec[7:8])));

  // loglikelihood
    target += reduce_sum(AR_MNLwSimPop_max, array_slice, splitsize,
                   sim_pop_int, sim_pop_price, morph_rowtocol,
                   share, 
                   log_price,  
                   b_dist, log_dist,
                   b_trend,
                   int_sku_over, int_sku_npl,
                   b_promo_exist, b_promo_npl, promo,
                   b_aware, log_aware,
                   b_morph_npl_self, b_morph_npl_hyp, b_morph_exist_hyp,
                   b_channel_exist, b_channel_npl, channel_bin,
                //   b_tier_sku, tier_price_sku,
                   per_base, per_lag, per_new, region,
                   skus_bin, lag_use,
                   wts, wts_att,
                   hill_g, hill_a
                  );

} // End Model