######################################
##  1. RUN: LOAD PACKAGES           ##
#####
library(readxl) # Read excel
library(cmdstanr) # Interface to Stan
library(posterior) # Process Stan draws
library(doParallel); library(parallel) # R Multi-threading for EB
source("https://raw.githubusercontent.com/klattery/Estimation-Ecosystem/master/EE_Functions_Stan2.R")
filename <- list() # Placeholder for filenames to import (data, covariates, excel specs)
set_cmdstan_path("/home/rstudio/cmdstan")
dir <- list(data = "/home/rstudio", # data files
            work = "/home/rstudio", # output and other work
            stanmodel = "/home/rstudio/StanCode", # code for Stan Model
            stanout = "/home/rstudio/StanCode" # Stan output 
)
#####

######### `MANDATORY INPUT` ###########################
##  *2. UPLOAD FILES to AWS and SPECIFY NAMES*       ##
##     VERIFY data_conjoint, specs_, data_cov        ##
#####
# Specify names for data file, Excel attribute/constraints specs, and output
out_prefix <- "MyOut" # YOUR name for output files (prefix)
filename$conjoint <- "MyData.csv/.rds" # YOUR .csv or .rds with extension
filename$specs <- "MyExcelSpecs.xlsx" # YOUR Coding and constraints file
filename$cov <- NULL # NULL if No Covariates. OPTIONAL file has id as 1st column

# 2) Read Data into R =============  
if (!file.exists(file.path(dir$data, filename$specs))){ # Read excel specs
  message(paste0("ERROR: CANNOT FIND YOUR Excel SPECS FILE: ", filename$specs))    
} else {
  specs_att_coding <- data.frame(read_xlsx(file.path(dir$data,filename$specs), sheet = "Att_Coding",
                                           col_types = c("text","text","numeric")))
  specs_pair_constraints <- data.frame(read_xlsx(file.path(dir$data,filename$specs), sheet = "Pair_Constraints",
                                                 col_types = c("text","numeric","numeric")))
  specs_cov_coding <- data.frame(read_xlsx(file.path(dir$data,filename$specs),
                                           sheet = "Cov_Coding", col_types = c("text","text")))
  message("\nREAD SPECS INTO R files:\n")
  message("specs_att_coding"); print(specs_att_coding)
  message("\nspecs_pair_constraints"); print(specs_pair_constraints)
  message("\nspecs_cov_coding"); print(specs_cov_coding);
}
data_conjoint <- read_csv_rds(dir$data, filename$conjoint, "as data_conjoint") # Reads in (csv or rds)
data_cov <- read_csv_rds(dir$data, filename$cov, "as data_cov") # Reads in (csv or rds)
#####

####################################################
##    3. RUN: CODE AND PREPARE DATA FOR STAN      ##
##       CHECK: PRINTED OUTPUT LOOKS RIGHT        ##
##       RECOMMEND: REVIEW CODE_MASTER_ANDCON.CSV ##
#####
col_id_task_dep <- c(1,2,ncol(data_conjoint)) # Columns for id, task, dep
# 3) Code and Prepare data_stan ============= 
indcode_spec <- indcode_spec_files(data_conjoint, specs_att_coding,
                                   specs_pair_constraints) # Code and constraints each attribute
# Here is where you would do any custom coding indcode_spec[[i]] <- "syntax"
indcode_list <- make_codefiles(indcode_spec) # Combine specifications above into one list
save_codemastercon(indcode_list,dir$work, out_prefix) # Save combined code_master and constraints 
data_stan <- prep_file_stan(idtaskdep = data_conjoint[,col_id_task_dep],
                            indcode_list)
if (!is.null(data_cov)){ # Code respondent covariates
  data_stan$i_cov <- code_covariates(data_cov, specs_cov_coding, data_stan$resp_id) 
  data_stan$P_cov <- ncol(data_stan$i_cov) # Num of coded parameters
}
#####

###################################################
##    4. RUN: (OPTIONALLY CHANGE DEFAULTS)       ##
#####
# 4) Modeling Options ============= 
# Specify multi-threading Stan and R
rm(indcode_spec); rm(indcode_list); gc()
threads_per_chain <- min(max(1,(detectCores() - 2)/2), round(.5 + data_stan$T/(1000)), 24)

#==================
# Modeling parameters. Defaults are usually fine
data_model <- list(
  iter_warmup = 400, # warmup of 400 is plenty
  iter_sampling = 400, # sampling of 400 is plenty
  df = 2,              # recommend df = 2 for Wishart
  prior_cov_scale = 1, # default 1
  splitsize = round(.5 + data_stan$T/(4 * threads_per_chain)), # Tasks per core
  agg_model = NULL, tag = NULL, ind = NULL
)
#####

###################################################
##   5. RUN: ESTIMATE MODEL                      ##
#####
# 5) Estimate Model ============= 
#####  Specify Stan Model 
stan_file <- "BaseHB_wPairCon_v2.1.stan" # Name of stan model in dir$stanmodel
stan_outname <- paste0(out_prefix, "_StanOut_", 
                       format(Sys.time(), '%Y%m%d-%H%M%S')) # Base Name of Stan Output files  

#####  Run Stan Model 
HB_model <- cmdstan_model(file.path(dir$stanmodel,stan_file), quiet = TRUE,
                          cpp_options = list(stan_threads = TRUE))
message_estimation(dir, stan_outname)
HB_fit <- HB_model$sample(modifyList(data_stan, data_model),
                          iter_warmup = data_model$iter_warmup,
                          iter_sampling = data_model$iter_sampling,
                          output_dir = dir$stanout,
                          output_basename = stan_outname, # set stan_outname above if changing
                          chains = 2,
                          parallel_chains = 2,
                          threads_per_chain = threads_per_chain,
                          save_warmup = TRUE,
                          refresh = 10,
                          adapt_delta = .8,
                          seed = 271,
                          init = .1,
                          show_messages = FALSE,
                          validate_csv = FALSE)
saveRDS(HB_fit, file.path(dir$work, "HB_fit.rds"))

#####  Check Convergence, Export Files  
if (max(HB_fit$return_codes() == 0)){
  checkconverge_export(stan_outname, dir$stanout, nchains = HB_fit$num_chains(),
                       vnames = colnames(data_stan$code_master), out_prefix, dir$work)
  process_utilities(data_stan, utilities, out_prefix, dir$work)
} else message("Stan Estimation Did not Finish")
#####

## OPTIONAL EMPIRICAL BAYES ##
#####
# 6) Optional Empirical Bayes =============
linux <- (Sys.info()[1] != "Windows")
if (linux) r_cores <- min(8,max(detectCores() -1,1)) # Number of cores to use for EB
if (!linux) r_cores <- min(4,max(detectCores() -1,1)) # Number of cores to use for EB
eb_betas_est(data_stan, draws_beta, colMeans(utilities), r_cores,
             out_prefix, dir$work, cov_scale = 1, linux)
#####
