#######################################################
#  1. Tools > Terminal > Terminal Options > Open Bash #
#  Run in Unix Terminal Ctrl-Alt-Enter                #
##################################################
# Multiplied attribute weights by .1

dir_data <- "C:/Users/K.Lattery/SKIM/Philip Morris - F6886 Advanced modelling update Nov 21/3. Modelling/01 RedWhite"
dir_out <- file.path(dir_data, "NoTrend")
dir_model <- "C:/Users/K.Lattery/Documents/GitHub/Stan/Test"
dir_draws <- "C:/StanRuns/PMI" # Where Stan stores draws.  Recommend a folder that does not sync

linux <- FALSE
if (linux) {
  dir_data <- gsub("C:/", "/mnt/c/", dir_data)
  dir_out <- gsub("C:/", "/mnt/c/", dir_out)
  dir_model <- gsub("C:/", "/mnt/c/", dir_model)
  dir_draws <- gsub("C:/", "/mnt/c/", dir_draws)
}

# Specify threads for your computer
threads = list(parallel_chains = 2,
               threads_per_chain = 15)

# Load data into R
data_list <- readRDS(file.path(dir_data, "data_list_AR_RW_wpromo.rds")) # Load data file
data_list <- readRDS(file.path(dir_data, "data_list_AR_RW_wpromo2021.rds")) # Load data file

bad_lag <- (data_list$per_lag <= 2) # Earliest lag period is 3
########################################
##    Special Red & White     ##########
bad_lag <- (data_list$per_lag <= 2) |
  (data_list$per_lag %in% (5:9)) |
  (data_list$per_new %in% (5:9)) 
#######################################

data_list$per_lag <- data_list$per_lag[!bad_lag]
data_list$per_new <- data_list$per_new[!bad_lag]
data_list$region <- data_list$region[!bad_lag]
data_list$T <- length(data_list$per_lag)
cbind(data_list$per_lag, data_list$per_new)

look <- t(data_list$aware[1,,])
look <- cbind(data_list$codes, look)
look <- data_list$code_master_data
#lag_dist_prev <- data_list$lag_dist
#new_dist_prev <- data_list$new_dist
#data_list$lag_dist <- log(exp(data_list$lag_dist) + 0) # Changed distribution
#data_list$new_dist <- log(exp(data_list$new_dist) + 0)
#plot(lag_dist_prev, data_list$lag_dist)

# Specify Constraints and Priors
P <- data_list$P
I <- 1000
dim(data_list$wts)
data_list$wts[1,dim(data_list$wts)[2],]
#for (i in 62:65){ # Create holdouts
#  data_list$wts[,i,] <- 0
#}
# Tuning weights for att model so later periods weighted much more
max_p <- dim(data_list$wts)[2]; min_p <-1
wts_period <- exp(2/(max_p - min_p) *(min_p:max_p - max_p))
wts_period <- wts_period/mean(wts_period)

########################################
##    Special Red & White     ##########
wts_period[5:9] <- 0
#######################################

#wts_period <- wts_period * .1
plot(1:length(wts_period), wts_period)
for (i in min_p:max_p){
  data_list$wts[,i,] <- wts_period[i]
  # data_list$wts_att[,i,] <- wts_period[i] * .25
  data_list$wts_att[,i,] <- wts_period[i] * 1
}
att_sigma <- rep(10, P-1)
as.matrix(colnames(data_list$sku_to_att))
att_sigma[46:55] <- 10 # Tier_format
att_sigma[56:61] <- 1 # Flavors
att_sigma[62:64] <- 2 # Colors
data.frame(att = colnames(data_list$sku_to_att),att_sigma)

sku_n  <- colSums(data_list$skus_bin[1,,]) +
  colSums(data_list$skus_bin[2,,]) +
  colSums(data_list$skus_bin[3,,]) +
  colSums(data_list$skus_bin[4,,]) +
  colSums(data_list$skus_bin[5,,])
sku_mean <- sku_n/5
dist_sigma <- (sku_mean - 1)/10
dist_sigma[dist_sigma > 1] <- 1
dist_sigma[dist_sigma < .05] <- .05

data_model <- list(
  per_base = max_p,
  P_int = data_list$P - 1,
  I = I,
  sim_pop_z = matrix(rnorm(I * P), nrow = P, ncol = I),
  price_sigma = 2,
  p_mu_mu = -1,
  p_mu_sigma = 2,
  trend_sigma = .1,
  dist_sigma = dist_sigma, # Was 1
  cov_diag_sigma = 2, # initial = 2
  att_sigma = att_sigma,
  ar_scale = 1,
  df = 2,
  cov_block = matrix(1, P, P),
  splitsize = round(.5 + data_list$T/(4 * threads[[2]])),
  lam_inv = 5, # For morphing max lam = .2
  data_tag = NULL,
  data_agg = NULL,
  code_master_data = NULL,
  codes_array = NULL
)

ls()
str(modifyList(data_list, data_model))
lapply(modifyList(data_list, data_model), function(x) sum(is.na(x)))
lapply(modifyList(data_list, data_model), function(x) sum(is.infinite(x)))

convert_3d <- function(x, logtrans = FALSE){
  dims <- dim(x)
  lev_sub <- lapply(1:dims[2], function(x) NA) # Create list with elements
  result <- lapply(1:dims[1], function(x) lev_sub)
  for (i in 1:dims[[1]]){ # fill in data_list
    for (j in 1:dims[[2]]){
      vec_one <- x[i,j,]
      if(logtrans){
        vec_one <- log(vec_one)
        vec_one[is.infinite(vec_one)] <- 0
      }
      result[[i]][[j]] <- vec_one
      
    }}
  return(result)
}
convert_3d_m <- function(x){
  dims <- dim(x)
  result <- lapply(1:dims[1], function(x) NA)
  for (i in 1:dims[[1]]){ # fill in data_list
    result[[i]] <- x[i,,]
  }
  return(result)
}

test_function <- function(junk){
  data_list$codes_array <- convert_3d(data_list$codes_array)
  data_list$price <- convert_3d(data_list$price, logtrans = TRUE)
  data_list$dist <- convert_3d(data_list$dist, logtrans = TRUE)
  data_list$aware <- convert_3d(data_list$aware)
  data_list$promo <- convert_3d(data_list$promo)
  data_list$share <- convert_3d(data_list$share)
  data_list$wts <- convert_3d(data_list$wts)
  data_list$wts_att <- convert_3d(data_list$wts_att)
  data_list$skus_bin <- convert_3d(data_list$skus_bin)
  data_list$morph_rowtocol <- convert_3d_m(data_list$morph_rowtocol) # Convert to array of matrices
}

library("cmdstanr") 
#HB_model <- cmdstan_model(file.path(dir_model, "DataFusion_PMI_AR_DiffMNP_ModBartDiagv2.stan"), quiet = TRUE, cpp_options = list(stan_threads = TRUE))
stan_file <- "Stan_Code_BlockCheckMax_Int10_popCC4_RW_newpriors2.stan" # Bchannel exist and npl
stan_file <- "Stan_Code_BlockCheckMax_Int10_popCC4_RW_newpriors2.stan" # Bchannel exist and npl

str(modifyList(data_list, data_model))

HB_model <- cmdstan_model(file.path(dir_model,stan_file), quiet = TRUE,
                          cpp_options = list(stan_threads = TRUE))
data_model$task_beg <- 1
data_model$task_end <- data_list$T
as.matrix(colnames(data_list$code_master_data))
data_model$nbrand <- 41 # Lst brand level
HB_MLE <- HB_model$optimize(modifyList(data_list, data_model), init = .2, seed = 0916,
                            refresh = 1, iter = 1000, output_dir =  dir_draws, threads = threads[[2]])
HB_MLE$output()
HB_model$print() # Optional to view

##############################################
# FOR OPTIMIZATION
# Get results
##############################################
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202202201654-1-1de050.csv"
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202202210946-1-2973d5.csv"
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202202211356-1-5b46ed.csv"
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202202211655-1-6bb3f2.csv"
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202203040907-1-43b736.csv"
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202203131640-1-3e1672.csv"
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202203131801-1-3dea9a.csv"
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202203141442-1-3508e9.csv" # Gaps 1,2
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202203150811-1-1cad07.csv" # Gaps 1,2,6
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202203160901-1-2035e9.csv" # Gaps 1,2,6 1 Master
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202203160901-1-2035e9.csv" # Gaps 1,2,6 1 Master
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202203171403-1-7a799d.csv"
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202203180855-1-052047.csv"
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202203181528-1-534824.csv" # MOd dist .1
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202203181743-1-691fa1.csv" # MOd dist .5
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202203191634-1-15b0f0.csv" # MOd dist .1 Prior Code 0
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202203192004-1-71e582.csv" # MOd dist .01 Prior Code 0
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202203192237-1-9351ff.csv" # MOd dist 0 (none) Prior Code 0
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202203200801-1-5f7902.csv" # dist 0, Prior Code 0, New wts
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart2.4-202203201020-1-17ddb6.csv" # dist 0, Prior Code 0, wts exp(5)
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart_promo-202203210839-1-55ded2.csv"
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart_promo-202203211410-1-4f062e.csv"
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart_promo-202203231300-1-9408aa.csv"
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart_promo-202203231801-1-806998.csv"
outname <- "Fusion_PMI_AR_Att_toSKU_corwdiag_bart_promo_foradj-202203240928-1-11e650.csv"
outname <- "Stan_Code_AWS_BlockCheckMax1-202205092133-1-988fb4.csv"
outname <- "Stan_Code_AWS_BlockCheckMax1.4_test-202205122146-1-2721d9.csv"
outname <- "Stan_Code_AWS_BlockCheckMax1.5-202205141350-1-059c84.csv"
outname <- "Stan_Code_AWS_BlockCheckMax1.6-202205150903-1-0d57a7.csv"
outname <- "Stan_Code_AWS_BlockCheckMax1.6-202205152028-1-81018e.csv"
outname <- "Stan_Code_AWS_BlockCheckMax1.7-202205161021-1-617cac.csv"
outname <- "Stan_Code_AWS_BlockCheckMax1.8-202205162329-1-207247.csv"
outname <- "Stan_Code_AWS_BlockCheckMax1.8-202205171429-1-6f9031.csv" # Smaller weights
outname <- "Stan_Code_AWS_BlockCheckMax1.9-202205172040-1-2391ff.csv"
outname <- "Stan_Code_AWS_BlockCheckMax1.10-202205181200-1-17a642.csv"
outname <- "Stan_Code_AWS_BlockCheckMax1.10-202205181952-1-17212c.csv"
outname <- "Stan_Code_AWS_BlockCheckMax1.10-202205192145-1-80852c.csv"
outname <- "Stan_Code_BlockCheckMax_Int-202205261531-1-6849ec.csv"
outname <- "Stan_Code_BlockCheckMax_IntAtt1-202205271056-1-163fc6.csv"
outname <- "Stan_Code_BlockCheckMax_Int3-202205281447-1-5adc11.csv"
outname <- "Stan_Code_BlockCheckMax_Int3-202205291209-1-830020.csv"
outname <- "Stan_Code_BlockCheckMax_Int4-202205301033-1-5dff18.csv"
outname <- "Stan_Code_BlockCheckMax_Int6_pop3_check-202206052148-1-72161d.csv"
outname <- "Stan_Code_BlockCheckMax_Int7_pop4-202206062005-1-0b010a.csv"
outname <- "Stan_Code_BlockCheckMax_Int8_pop-202206090804-1-2f6709.csv"
outname <- "Stan_Code_BlockCheckMax_Int10_pop-202206111002-1-5e1201.csv"
outname <- "Stan_Code_BlockCheckMax_Int10_popCC_test-202206241127-1-936128.csv" # normal(0,2)
outname <- "Stan_Code_BlockCheckMax_Int10_popCC2-202206242019-1-2b8961.csv" # normal(0,5)
outname <- "Stan_Code_BlockCheckMax_Int10_popCC2-202206250911-1-2d1f92.csv" # normal(-.1, 1)
outname <- "Stan_Code_BlockCheckMax_Int10_popCC2-202206252019-1-11f1cd.csv" # normal(-.2, 1)
outname <- "Stan_Code_BlockCheckMax_Int10_popCC2-202206260826-1-47734c.csv" # normal(0,5)
outname <- "Stan_Code_BlockCheckMax_Int10_popCC3-202206272035-1-6f2942.csv"
outname <- "Stan_Code_BlockCheckMax_Int10_popCC4-202206291455-1-73c756.csv"
outname <- "Stan_Code_BlockCheckMax_Int10_popCC4_collapse-202206301004-1-0da812.csv"
outname <- "Stan_Code_BlockCheckMax_Int10_popCC4_collapse_Bristol_test-202207142052-1-3599da.csv"
outname <- "Stan_Code_BlockCheckMax_Int10_popCC4_collapse_Bristol_3per-202207230935-1-746525.csv"
outname <- "Stan_Code_BlockCheckMax_Int10_popCC4_collapse_Bristol_3per_newpriors-202207310858-1-8c170a.csv"
outname <- "Stan_Code_BlockCheckMax_Int10_popCC4_collapse_Bristol_3per_newpriors2-202208032107-1-7554dc.csv"
outname <- "Stan_Code_BlockCheckMax_Int10_popCC4_collapse_Bristol_3per_newpriors2-202208041550-1-7aa3a5.csv"
outname <- "Stan_Code_BlockCheckMax_Int10_popCC4_collapse_Bristol_3per_newpriors3-202208042228-1-2056ce.csv"
outname <- "Stan_Code_BlockCheckMax_Int10_popCC4_collapse_Bristol_3per_newpriors3a-202208061849-1-7904ab.csv"
outname <- "Stan_Code_BlockCheckMax_Int10_popCC4_collapse_Bristol_3per_newpriors4-202208111045-1-758b9f.csv"
outname <- "Stan_Code_BlockCheckMax_Int10_popCC4_RW_newpriors2-202208171020-1-0d7bf9.csv"
outname <- "Stan_Code_BlockCheckMax_Int10_popCC4_RW_newpriors2-202208192119-1-3c236e.csv"
outname <- "Stan_Code_BlockCheckMax_Int10_popCC4_RW_newpriors2-202208241656-1-3b538d.csv"
outname <- "Stan_Code_BlockCheckMax_Int10_popCC4_RW_newpriors2-202209021154-1-90bb8c.csv"
outname <- "Stan_Code_BlockCheckMax_Int10_popCC4_RW_newpriors2-202209021804-1-4e1a14.csv"
outname <- "Stan_Code_BlockCheckMax_Int10_popCC4_RW_newpriors2-202209051030-1-2a2a39.csv"

#b_price <- as.vector(read_cmdstan_csv(file.path(dir_draws,outname), variables = "b_price")$point_estimates)
cov_chol_v <- as.vector(read_cmdstan_csv(file.path(dir_draws,outname), variables = "cov_chol")$point_estimates)
b_dist <- as.vector(read_cmdstan_csv(file.path(dir_draws,outname), variables = "b_dist")$point_estimates)
b_trend <- as.vector(read_cmdstan_csv(file.path(dir_draws,outname), variables = "b_trend")$point_estimates)
p_mu <- as.vector(read_cmdstan_csv(file.path(dir_draws,outname), variables = "p_mu")$point_estimates)
sd_diag <- as.vector(read_cmdstan_csv(file.path(dir_draws,outname), variables = "sd_diag")$point_estimates)
b_promo <- as.vector(read_cmdstan_csv(file.path(dir_draws,outname), variables = "b_promo")$point_estimates)
b_morph_npl_self <- as.vector(read_cmdstan_csv(file.path(dir_draws,outname), variables = "b_morph_npl_self")$point_estimates)
b_morph_npl_hyp <- as.vector(read_cmdstan_csv(file.path(dir_draws,outname), variables = "b_morph_npl_hyp")$point_estimates)
b_morph_exist_hyp <- as.vector(read_cmdstan_csv(file.path(dir_draws,outname), variables = "b_morph_exist_hyp")$point_estimates)

b_aware <- as.vector(read_cmdstan_csv(file.path(dir_draws,outname), variables = "b_aware")$point_estimates)
b_attributes_over <- as.vector(read_cmdstan_csv(file.path(dir_draws,outname), variables = "b_attributes_over")$point_estimates)
b_attributes_npl <- as.vector(read_cmdstan_csv(file.path(dir_draws,outname), variables = "b_attributes_npl")$point_estimates)
b_channel_exist <- as.vector(read_cmdstan_csv(file.path(dir_draws,outname), variables = "b_channel_exist")$point_estimates)
b_channel_npl <- as.vector(read_cmdstan_csv(file.path(dir_draws,outname), variables = "b_channel_npl")$point_estimates)

look <- data.frame(colnames(data_list$prior_cor), round(c(b_attributes_over,p_mu),5), round(c(b_attributes_npl, p_mu),5))
exp(b_promo)
look <- data.frame(colnames(data_list$prior_cor), c(b_attributes_over, p_mu))
look <- data.frame(colnames(data_list$prior_cor), c(colnames(data_list$sku_to_att),"logprice"))
check <- data.frame(data_list$codes, data_list$sku_to_att)

# Look at covariance, correlation
cov_chol <- matrix(cov_chol_v, data_list$P)
#cov_m <- diag(sd_diag) %*% cov_chol %*% t(cov_chol) %*% diag(sd_diag)
cov_m <- cov_chol %*% t(cov_chol)

klower <- lower.tri(cov_m)
plot(as.vector((data_list$prior_cor)[klower]), as.vector((cov_m)[klower]))
plot(as.vector(cov2cor(data_list$prior_cor)[klower]), as.vector(cov2cor(cov_m)[klower]), xlab = "Prior Correlation", ylab = "Correlation from Estimated Real World Covariance", col = "blue")

#b_dist <- rep(1,P)
round(b_trend,3)
round(b_dist,3)
#hist(b_price, breaks = 30)

len_short <- data_model$P_int
vnames <- colnames(data_list$prior_cor)
vnames_short <- vnames[1:len_short]

sim_pop_est <- t(cov_chol %*% data_model$sim_pop_z)
plot(as.vector(cov_m), as.vector(cov(sim_pop_est)))
cor_est <- cor(sim_pop_est)
klower <- lower.tri(data_list$prior_cor)
plot(as.vector(cov2cor(data_list$prior_cor)[klower]), as.vector(cor_est[klower]))
hist(as.vector(cov2cor(data_list$prior_cor)[klower]), breaks = 30)
hist(as.vector(cor_est[klower]), breaks = 30)

z_pop <- matrix(rnorm(10000 * ncol(cov_chol)), nrow = ncol(cov_chol))
sim_pop <- t(cov_chol %*% z_pop)
plot(as.vector(cov2cor(data_list$prior_cor)[klower]), as.vector(cor(sim_pop)[klower]))

sim_pop_price_raw <- sim_pop[,ncol(sim_pop)]
sim_pop_price <- -.2 * log(1 + exp(1/(-.2) * (p_mu + sim_pop_price_raw)))
mean(sim_pop_price); sd(sim_pop_price); sd(sim_pop_price_raw)
hist(sim_pop_price, breaks = 100)
plot(sim_pop_price_raw, sim_pop_price)

sim_pop_int <- sim_pop[,1:len_short]
sim_pop_sku <- sim_pop_int %*% t(data_list$sku_to_att) # create sku covariance for testing
colnames(sim_pop_int) <- vnames_short
colnames(sim_pop_sku) <- data_list$skus_att$variant
round(colMeans(sim_pop_sku),3); hist(colMeans(sim_pop_sku))
#sim_pop_sku <- scale(sim_pop_sku, scale = FALSE)
apply(sim_pop_sku, 2, sd) 
sim_pop_full <- sim_pop %*% t(data_list$code_master)

names(b_attributes_over) <- vnames_short
names(b_attributes_npl) <- vnames_short
code_master_xprice <- data_list$code_master_data
code_master_xprice <- code_master_xprice[1:nrow(code_master_xprice)-1,1:len_short]
b_attributes_npl_full <- as.vector(b_attributes_npl %*% t(code_master_xprice))
b_attributes_over_full <- as.vector(b_attributes_over %*% t(code_master_xprice))
sku_to_att_full <- data_list$sku_to_att %*% t(code_master_xprice)
names(b_attributes_npl_full) <- rownames(code_master_xprice)
names(b_attributes_over_full) <- rownames(code_master_xprice)
b_attributes_compare <- as.matrix(cbind(b_attributes_over_full, b_attributes_npl_full))
write.csv(b_attributes_compare, file = file.path(dir_out, "b_attributes_compare.csv"))

result_optim_att_sku <- list(sku_to_att = data_list$sku_to_att,
                             sku_to_attwcode = cbind(code_variant = data_list$codes, data_list$sku_to_att),
                             sku_to_attwcode_full = cbind(code_variant = data_list$codes, data_list$sku_to_att_full),
                             sku_betas_wcode = data.frame(code_variant = data_list$codes,
                                                          b_dist = round(b_dist,5),
                                                          b_trend = round(b_trend,5)),
                             b_attributes_npl = b_attributes_npl, b_attributes_npl_full = b_attributes_npl_full,
                             b_attributes_over = b_attributes_over, b_attributes_over_full = b_attributes_over_full,
                             sim_pop_int = round(sim_pop_int,5), sim_pop_price = round(sim_pop_price,5), sim_pop_full = round(sim_pop_full,5), sim_pop_sku = sim_pop_sku,
                             code_master = data_list$code_master, cov_chol = cov_chol, b_promo = b_promo, b_aware = b_aware,
                             b_morph_npl_self = b_morph_npl_self, b_morph_npl_hyp = b_morph_npl_hyp, b_morph_exist_hyp =  b_morph_exist_hyp,
                             b_channel_exist = b_channel_exist, b_channel_npl = b_channel_npl,
                             ar_scale =1,
                             per_base = data_model$per_base)
betas_new <- (colMeans(result_optim_att_sku$sku_betas_wcode[,-1]) + c(1,0))/2
result_optim_att_sku$sku_betas_new <- betas_new

saveRDS(result_optim_att_sku, file = file.path(dir_out, "result_optim_att_sku_RW_3per_newprior4_2021.rds"))
