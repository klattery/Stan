# Uses model with price in population
library(readxl)
library(sqldf)
library(shiny)

dir_data <- "C:/Users/K.Lattery/SKIM/Philip Morris - F6886 Advanced modelling update Nov 21/3. Modelling/01 GT IMS top3"
dir_work <- file.path(dir_data, "Model3/Coding1_wInt8_PricePop")
dir_base <- "C:/Users/K.Lattery/SKIM/Philip Morris - 3. Prework/3 final data"
data_tag <- file.path(dir_data, "IMS_Top3_PanRussia_Mar2022_wpromo.rds") # Pan Russia

data_full <- readRDS(data_tag) # Stacked data
check <- unique(data_full[,c(1,2)])
data_use <- data_full[data_full$dist> .05,]
codes_data <- sort(unique(data_use$code_variant))
model_pars <- readRDS(file.path(dir_work,"result_optim_att_sku_hold_promo65_wint.rds"))
#model_pars$per_base <- 60 # I forget to store it in result

#############################################################
#1) Match our skus to those matching attributes in our model
############################################################
att_decompose <- readRDS(file.path(dir_base,"02 Modelling/EvokedwAtts","sku_att_all3.rds"))
kmatch <- match(codes_data, att_decompose$variant)
sum(is.na(kmatch)) # Check sum = 0, all our skus are in database
if (sum(is.na(kmatch)) > 0){
  message("Some skus not in database of attributes")
  skus_either[is.na(kmatch)]
} 
skus_att <- att_decompose[kmatch[!is.na(kmatch)],] # These are skus that match our attribute database 

# Match to attributes in our model
kmatch <- colnames(skus_att) %in% c(colnames(skus_att)[1:2], colnames(model_pars$sku_to_att))
badcols <- skus_att[,!kmatch,drop = FALSE]
colnames(badcols)
if (ncol(badcols)>0){
  badrows <- apply(badcols, 1, function(x){
    result <- TRUE
    if (min(x) == 0 & max(x) == 0) result <- FALSE
    return (result)
  })
  message("Dropping the following skus because they have attributes not in model:")
  print(skus_att[badrows,c(1,2)])
  skus_att <- skus_att[!badrows,kmatch]
}
sum(colnames(skus_att[,c(-1:-2)]) != colnames(model_pars$sim_pop))
#cbind(colnames(skus_att[,c(-1:-2)]), colnames(model_pars$sim_pop))

data_use <- data_use[data_use$code_variant %in% skus_att$variant,]
codes_use <- sort(unique(data_use$code_variant))
P_sku <- length(codes_use)

###################################################
#2) Convert data to block form region x period x sku code
###################################################
regions <- sort(unique(data_use$region_num)); N_regions <- length(regions)
periods <- sort(unique(data_use$period)); N_periods <- length(periods)
regions; periods # SHOULD BE 1:N
codes_array <- array(data = 0, dim = c(N_regions, length(periods), P_sku))
var_initial <- c("codes_array", "price", "dist", "aware", "promo", "share", "wts", "wts_att", "skus_bin")
data_list <- lapply(var_initial, function(x) assign(x, codes_array)) # Create list with elements
names(data_list) = var_initial
data_list$aware[] <- 0
for (i in regions){ # fill in data_list
  for (j in periods){
    x <- data_use[data_use$region_num == i & data_use$period == j ,]
    kmatch <- match(x$code_variant, codes_use)
    data_list$codes_array[i,j,kmatch] = x$code_variant # diagnostic check only
    data_list$price[i,j,kmatch] = x$prices_use
    data_list$dist[i,j,kmatch]  = x$dist
    data_list$aware[i,j,kmatch] = x$aware
    data_list$promo[i,j,kmatch] = x$promo
    data_list$share[i,j,kmatch] = x$packs/sum(x$packs)
    data_list$wts[i,j,kmatch]   = 1
    data_list$wts_att[i,j,kmatch] = 1
    data_list$skus_bin[i,j,kmatch] = 1
  }
}
data_list$N_regions <- N_regions
data_list$N_periods <- N_periods
data_list$P_sku <- P_sku
data_list$codes <- codes_use
data_list <- c(data_tag = data_tag, data_list)
data_list$sku_to_att <- as.matrix(skus_att[,c(-1,-2)])

###################################################
# 3. Add Morphing
##################################################
dir_morph <- "C:/Users/K.Lattery/SKIM/Philip Morris - F6886 Advanced modelling update Nov 21/2. Data/2 DP SKIM/1 All channels merged"
data_morph <- data.frame(read_xlsx(file.path(dir_morph, "All_channels_2021_total Russia.xlsx"),
                                   range = "info!P3:V28"))
morph_sku1 <- match(data_morph$variant1, data_list$codes)
morph_sku2 <- match(data_morph$variant2, data_list$codes)
bad <- is.na(morph_sku1) | is.na(morph_sku2)
data_morph <- cbind(data_morph, morph_sku1, morph_sku2)[!bad,]

# For each period construct sku x sku matrix
get_morph_matrix <- function(pnum, data_morph, codes){
  nskus <- length(codes)
  result <- matrix(0, nskus, nskus)
  morph_period <- data_morph[pnum >= data_morph$period & pnum <= data_morph$period + 12,] # Morph started and less than 1 year
  if (nrow(morph_period) > 0){
    sku_from <- match(morph_period$variant1, codes)
    sku_to <- match(morph_period$variant2, codes)
    result[cbind(sku_from, sku_to)] <- 1
  }
  return(result)
}

morph_to_list <- lapply(1:data_list$N_periods, function(i) get_morph_matrix(i, data_morph, data_list$codes))
morph_rowtocol <- array(dim = c(length(morph_to_list),dim(morph_to_list[[1]]))) # convert list to 3d Matrix
for (i in 1:length(morph_to_list)){
  morph_rowtocol[i,,] <- morph_to_list[[i]]
}
data_list$morph_rowtocol <- morph_rowtocol
data_list$data_agg <- data_use # Keeping for reference
# saveRDS(data_list, file = file.path(dir_work, "data_list_AR_IMS_Top3_wpromo65_block.rds")) 

###########################################################
# 4. Fill in sku_betas for new products not in estimation
##########################################################
# kmatch <- match(model_pars$sku_betas_wcode$code_variant, data_list$codes)
# price_reg <- data.frame(cbind(b_price = model_pars$sku_betas_wcode$b_price, model_pars$sku_to_att))
# fit_linear <- lm(b_price~., data = price_reg)
# summary(fit_linear)
# #plot(price_reg$b_price, fit_linear$fitted.values)
# p_betas <- fit_linear$coefficients
# p_betas[is.na(p_betas)] <- 0
# pred_p <- as.matrix(cbind(1, price_reg[,-1])) %*% p_betas
# plot(price_reg$b_price, pred_p)
# hist(price_reg$b_price - pred_p)

# Add new sku betas 
new_skus <- setdiff(data_list$codes, model_pars$sku_betas_wcode$code_variant) # new skus without betas
kmatch <- match(new_skus, data_list$codes)
sku_betas_add <- cbind(code_variant = new_skus, b_dist = model_pars$sku_betas_new[1], b_trend = 0)
sku_betas_wnew <- rbind(model_pars$sku_betas_wcode, sku_betas_add)
kmatch <- match(data_list$codes, sku_betas_wnew$code_variant)
sku_betas_wnew <- sku_betas_wnew[kmatch,]
if(sum(sku_betas_wnew$code_variant != data_list$codes) > 0){
  message("FATAL ERROR!! SKU betas do not match skus in data.  Model will not forecast.")
} # Check 0
data_list$sku_betas_match <- sku_betas_wnew

###################################################
# 5. Forecasts
##################################################
maxexp_u_whatif<- function(sim_pop_one, mu_self, mu_whatif, morph_rowtocol, morph_receive, mult_self_b, mult_from_b, P_sku){
  # morph_receive is just col_max(morph_rowtocol) item receives morphing
  exp_util <- mult_self_b * exp(mu_self + sim_pop_one) # vector mult
  exp_util_whatif <- exp(mu_whatif + sim_pop_one)
  for (i in 1:P_sku){
    if (morph_receive[i] >= .99){
      # max_morph <- max(morph_rowtocol[,i] * exp_util_from) # Max from items morphing into i
      exp_util[i] <- max(exp_util[i],
                         mult_from_b[i] * max(morph_rowtocol[,i] * exp_util_whatif))
    }
  }
  return(exp_util)
}
# maxexp_u_whatif_one <- function(sim_pop_one, mu_self, mu_whatif, morph_rowtocol, morph_receive, mult_self_b, mult_from_b, P_sku){
#   # morph_receive is just col_max(morph_rowtocol) item receives morphing
#   exp_util <- mult_self_b * exp(mu_self + sim_pop_one) # vector mult
#   exp_util_whatif <- exp(mu_whatif + sim_pop_one)
#   for (i in 1:P_sku){
#     if (morph_receive[i] >= .99){
#       # max_morph <- max(morph_rowtocol[,i] * exp_util_from) # Max from items morphing into i
#       exp_util[i] <- max(exp_util[i],
#               mult_from_b[i] * max(morph_rowtocol[,i] * exp_util_whatif))
#     }
#   }
#   return(cbind(exp_util, exp_util_whatif, skus_new * exp_util_final))
# }


log_w0 <- function(index_binary, x){
  return(index_binary * log(x + 1.0000000001 - index_binary));
  #index = 1 returns log(x), index = 0 returns 0
}
col_max <- function(x){
  result <- apply(x, 2, max)
  return(result)
}
row_max <- function(x){
  result <- apply(x, 1, max)
  return(result)
}
# data transforms
log_price <- data_list$price
log_dist <- data_list$dist
log_aware <- data_list$aware
for (i in 1:data_list$N_regions){
  for (j in 1:data_list$N_periods){
    log_price[i,j,] <- log_w0(data_list$skus_bin[i,j,], data_list$price[i,j,] * .01); 
    log_dist[i,j,] <- log_w0(data_list$skus_bin[i,j,], data_list$dist[i,j,] * 10);
    log_aware[i,j,] <-  log(data_list$aware[i, j,] + .05) - log(1.05);
  }
}
morph_rowtocol <- list()
for (j in 1:data_list$N_periods){
  morph_rowtocol[[j]] <- data_list$morph_rowtocol[j,,]
}

per_new1 <- c(4:60) # Prior stuff
per_new2 <- 61:72 # Forecasts
per_lag1 <- per_new1 - 3
per_lag2 <- rep(60, length(per_new2))
per_new <- c(per_new1, per_new2)
per_lag <- c(per_lag1, per_lag2)
T <- length(per_new)
region <- rep(1, T)
cbind(per_lag, per_new, region)

AR_MNLwSimPop_max <- function(t,
                              sim_pop_int, sim_pop_price, morph_rowtocol,
                              share, 
                              b_price, log_price,
                              b_dist,  log_dist,
                              b_trend,
                              int_sku_over, int_sku_npl,
                              b_promo, promo,
                              b_aware, log_aware,
                              b_morph_npl_self, b_morph_npl_hyp, b_morph_exist_hyp,
                              per_base, per_lag, per_new,  region,
                              skus_bin){
  npop <- ncol(sim_pop_int)
  P_sku <- nrow(sim_pop_int)
  p_lag <- per_lag[t]; #lag period for task
  p_new <- per_new[t]; # new/forecast period for task
  region_t <- region[t];
  skus_new <- skus_bin[region_t, p_new,]; # // skus in new/forecast period: absolute
  skus_lag <- skus_bin[region_t, p_lag,]; # // skus in lag period: absolute
  skus_over <- skus_new * skus_lag; # // skus in both periods: relative
  skus_npl <- skus_new - skus_over; #// skus that are new product launches vs lag: relative
  mu_int <- (skus_over * int_sku_over) + (skus_npl * int_sku_npl)
  u_global_lag <- mu_int + b_trend * (p_lag - per_base) +
    b_dist  * log_dist[region_t, p_lag,]  +
    b_promo  * promo[region_t, p_lag,] +
    b_aware *  log_aware[region_t, p_lag,];
  
  u_global_fore <- mu_int + b_trend * (p_new - per_base) +
    b_dist  * log_dist[region_t, p_new,]  +
    b_promo  * promo[region_t, p_new,] +
    b_aware *  log_aware[region_t, p_new,];
  
  share_lag_over <- (share[region_t, p_lag,] * skus_over)/sum(share[region_t, p_lag,] * skus_over);
  #    check <- cbind(data_list$codes, skus_over, skus_npl, share_lag_over, u_global_lag, u_global_fore)
  #///////////////////
  #  // Simulations
  #////////////////////
  #  // Declare some storage
  task_prob_sum_new <- rep(0,P_sku);
  task_prob_sum_lag <- rep(0,P_sku);
  #sim_share_fore <- rep(0,P_sku);
  #sim_share_lag <- rep(0,P_sku);
  #vector[P_sku] U_final;
  #vector[P_sku] pred_new;
  #vector[P_sku] exp_util; // placeholder for e^U calcs
  
  #// Forecast
  #morph_rowtocol_clean <- diag_post_multiply(morph_rowtocol[p_new,,], skus_new); #//  Remove any skus that no longer exist
  morph_rowtocol_clean <- morph_rowtocol[p_new,,] %*% diag(skus_new); #//  Remove any skus that no longer exist
  #// util of sku @ price, dist, etc of sku morphing into
  u_global_whatif <-  int_sku_over +
    morph_rowtocol_clean %*% (b_trend * (p_new - per_base)) +
    morph_rowtocol_clean %*% (b_dist  * log_dist[region_t, p_new,])
  
  mu_whatif <- row_max(morph_rowtocol_clean) * u_global_whatif;     
  
  morph_receive <- col_max(morph_rowtocol_clean); #// column sums 0/1
  morph_npl <- morph_receive * skus_npl;    
  morph_exist <- morph_receive * (1-skus_npl);
  mult_self_b <- (skus_new - morph_npl) + morph_npl * b_morph_npl_self;
  mult_from_b <- morph_npl * b_morph_npl_hyp + morph_exist * b_morph_exist_hyp;
  
  log_price_t <- log_price[region_t, p_new,]
  log_price_t_whatif <- morph_rowtocol_clean %*% log_price_t;
  # check <- cbind(data_list$codes, skus_over, skus_npl, mult_self_b, mult_from_b, u_global_fore, mu_whatif)
  for (i in 1:npop){
    exp_util <- skus_new * maxexp_u_whatif(
      sim_pop_int[,i],
      u_global_fore + (sim_pop_price[i] * log_price_t), # change from v4
      mu_whatif + (sim_pop_price[i] * log_price_t_whatif), # change from v4
      morph_rowtocol_clean, morph_receive,
      mult_self_b, mult_from_b, P_sku);
    task_prob_sum_new <- task_prob_sum_new + exp_util/sum(exp_util); #// key: sum of shares over pop            
  }
  # sum(task_prob_sum_new)
  # i <- 4
  # look <- cbind(data_list$codes, skus_new, maxexp_u_whatif_one(
  #   sim_pop_int[,i] + (sim_pop_price[i] * log_price_t),
  #   u_global_fore, mu_whatif,
  #   morph_rowtocol_clean, morph_receive,
  #   mult_self_b, mult_from_b, P_sku));
  
  # Repeat for Lag  
  # Remove any skus not in overlap & that do not contribute to *FORECAST* (morph_self_clean forecast)
  #morph_rowtocol_clean <- diag_post_multiply(morph_rowtocol[p_lag,,], skus_over * morph_receive); #// yes this refers to forecast
  morph_rowtocol_clean <- morph_rowtocol[p_lag,,] %*% diag(skus_over * morph_receive); #// yes this refers to forecast
  # util of sku @ price, dist, etc of sku morphing into
  u_global_whatif =  int_sku_over +
    morph_rowtocol_clean %*% (b_trend * (p_lag - per_base)) +
    morph_rowtocol_clean %*% (b_dist  * log_dist[region_t, p_lag,])
  mu_whatif <- row_max(morph_rowtocol_clean) * u_global_whatif; 
  
  morph_receive <- col_max(morph_rowtocol_clean); #// column sums 0/1
  morph_npl <- morph_receive * skus_npl;    
  morph_exist <- morph_receive * (1-skus_npl);
  mult_self_b <- (skus_new - morph_npl) + morph_npl * b_morph_npl_self;
  mult_from_b <- morph_npl * b_morph_npl_hyp + morph_exist * b_morph_exist_hyp;
  
  log_price_t <- log_price[region_t, p_new,]
  log_price_t_whatif <- morph_rowtocol_clean %*% log_price_t;
  for (i in 1:npop){
    exp_util = skus_over * maxexp_u_whatif(
      sim_pop_int[,i] + (sim_pop_price[i] * log_price_t),
      u_global_lag, mu_whatif,
      morph_rowtocol_clean, morph_receive,
      mult_self_b, mult_from_b, P_sku);
    task_prob_sum_lag = task_prob_sum_lag + exp_util/sum(exp_util); #// key: sum of shares over pop 
  }
  
  # Compute final Utility and LogLikelihood
  sim_share_fore <- task_prob_sum_new/npop;
  sim_share_lag <- task_prob_sum_lag/npop;
  U_final1 <- skus_over * share_lag_over * sim_share_fore / (sim_share_lag + 1 - skus_over)
  U_final2 <- (1-skus_over) * sim_share_fore
  lag_ratio <- share_lag_over/(sim_share_lag + 1 - skus_over) # diagnostic want ~1
  U_final <- skus_over * share_lag_over * sim_share_fore / (sim_share_lag + 1 - skus_over) +
    (1-skus_over) * sim_share_fore;
  pred_new <- U_final/sum(U_final);
  
  result <- data.frame(data_list$codes, t = rep(t, P_sku), p_lag = rep(p_lag, P_sku), p_new = rep(p_new,P_sku), skus_over, skus_new,
                       lag_share = share[region_t, p_lag,], obs_share = share[region_t, p_new,],
                       pred_new, sim_share_lag, sim_share_fore, share_lag_over, lag_ratio)
  return(result);
}

forecast_all <- do.call(rbind, lapply(1:T, function(t){
  result <- AR_MNLwSimPop_max(
    t = t,
    sim_pop_int = data_list$sku_to_att %*% t(model_pars$sim_pop_int),
    sim_pop_price = model_pars$sim_pop_price,
    morph_rowtocol = data_list$morph_rowtocol,
    share = data_list$share,
    log_price = log_price,
    b_dist = data_list$sku_betas_match$b_dist,
    log_dist = log_dist,
    b_trend = data_list$sku_betas_match$b_trend,
    int_sku_over = data_list$sku_to_att %*% model_pars$b_attributes_over,
    int_sku_npl = data_list$sku_to_att %*% model_pars$b_attributes_npl,
    b_promo = model_pars$b_promo,
    promo = data_list$promo,
    b_aware =  model_pars$b_aware,
    log_aware = log_aware,
    b_morph_npl_self = model_pars$b_morph_npl_self,
    b_morph_npl_hyp = model_pars$b_morph_npl_hyp,
    b_morph_exist_hyp = model_pars$b_morph_npl_hyp,
    per_base = model_pars$per_base,
    per_lag = per_lag,
    per_new = per_new,
    region = region,
    skus_bin = data_list$skus_bin
  )
  result <- data.frame(codes = data_list$codes, result,
                       price_lag = data_list$price[region[t],per_lag[t],],
                       price_new = data_list$price[region[t],per_new[t],],
                       dist_lag = data_list$dist[region[t],per_lag[t],],
                       dist_new = data_list$dist[region[t],per_new[t],]
  )
  
  return(result)  
}))
forecast_new <- forecast_all[forecast_all$skus_new == 1,]
plot(forecast_new$obs_share, forecast_new$sim_share_fore, col = "orange") # sim shares
plot(forecast_new$obs_share, forecast_new$pred_new, col = "blue")
kfilter <- (forecast_new$skus_over == 0)
points(forecast_new$obs_share[kfilter], forecast_new$pred_new[kfilter], col = "red") # NPL
write.csv(forecast_new, file = file.path(dir_work, "forecast_new.csv"), row.names = FALSE)


#####################################
## Debug single task
t <- 68 # Nov 2021
sim_pop_int <- data_list$sku_to_att %*% t(model_pars$sim_pop_int)
sim_pop_price <- model_pars$sim_pop_price
morph_rowtocol <- data_list$morph_rowtocol
share <- data_list$share
log_price <- log_price
b_dist <- data_list$sku_betas_match$b_dist
log_dist <- log_dist
b_trend <- data_list$sku_betas_match$b_trend
int_sku_over <- data_list$sku_to_att %*% model_pars$b_attributes_over
int_sku_npl <- data_list$sku_to_att %*% model_pars$b_attributes_npl
b_promo <- model_pars$b_promo
promo <- data_list$promo
b_aware <-  model_pars$b_aware
log_aware <- log_aware
b_morph1 <- model_pars$b_morph1
b_morph2 <- model_pars$b_morph2
per_base <- model_pars$per_base
per_lag <- per_lag
per_new <- per_new
region <- region
skus_bin <- data_list$skus_bin
