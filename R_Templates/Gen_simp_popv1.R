gen_sim_pop <- function(alpha, cov, I, draw_alpha = FALSE, draw_cov = TRUE, con_sign = NULL){
  # Simulate population
  P <- length(alpha)
  if (P != ncol(cov)){
    stop("Fatal Error: alpha must have same length as columns of cov")
  }
  if (draw_alpha) alpha <- alpha + rnorm(length(alpha), 0 , 1)
  if (draw_cov) cov <- rWishart(1,P+2,cov/(P + 2))[,,1]
  cov_chol <- t(chol(cov)) # Lower triangular
  z <- matrix(rnorm(P * I), nrow = P, ncol = I)
  sim_pop <- sweep(t(cov_chol %*% z), 2, alpha, "+")
  if (!is.null(con_sign)){
    if (P != length(con_sign)){
      stop("Fatal Error: con_sign must have same length as columns of cov, alpha")
    }
    for (j in (1:P)[abs(con_sign)> 0]){
      x <- sim_pop[,j]
      mysign <- sign(con_sign[j])
      bad_x <- sign(x) != mysign
      sim_pop[bad_x,j] <- mysign * runif(sum(bad_x), 0, .3);
    }
  }
  return(sim_pop)
}  


sim_MNL <- function(data_stan, utilities){
  result <- do.call(rbind, lapply(1:data_stan$T,
            function(t){
              U <- (data_stan$ind[data_stan$start[t]:data_stan$end[t],] %*%
                      utilities[data_stan$task_individual[t],])
              pred <- exp(U)/sum(exp(U))
              choice <- maxU <- rep(0, length(U))
              choice[which.max(U + -log(-log(runif(length(U)))))] <- 1
              maxU[which.max(U)] <- 1
              return(data.frame(U= U, pred = pred, choice_noerr = maxU, choice_werr = choice))
            }
            ))
  return(result)
} 

sim_utility <- sim_pop[1,]
est_utility <- utilities[1,]
id <- 1
hold_tasks <- c(125:150, 250, 300, 800)
MNL_holdout_compare_one <- function(data_stan, sim_utility, est_utility, id, hold_tasks){
  # Compare simulated and estimated utility (1 respondent) on holdout tasks
  result <- do.call(rbind, lapply(hold_tasks,
                                  function(t){
                                    X <- data_stan$ind[data_stan$start[t]:data_stan$end[t],]
                                    U <- X %*% sim_utility
                                    sim_choice <- est_choice <- rep(0, length(U))
                                    sim_pred <- exp(U)/sum(exp(U))
                                    sim_choice[which.max(U)] <- 1
                                    U <- X %*% est_utility
                                    est_pred <- exp(U)/sum(exp(U))
                                    est_choice[which.max(U)] <- 1                                   
                                    return(data.frame(id = id, task = t,
                                                      sim_pred = sim_pred, sim_choice = sim_choice,
                                                      est_pred = est_pred, est_choice = est_choice))
                                  }
  ))
  return(result)
} 

utilities <- sim_pop + matrix(rnorm(600 * 21)/5, 600, 21)
plot(as.vector(sim_pop), as.vector(utilities), xlab = "Simulated Population Utilities", ylab = "Estimated Utilities")
plot(as.vector(cor(sim_pop)), as.vector(cor(utilities)), xlab = "Simulated Population Cor", ylab = "Estimated Correlation")
plot(as.vector(cov(sim_pop)), as.vector(cov(utilities)), xlab = "Simulated Population Cov", ylab = "Estimated Covaraince")

# Get holdouts each respondent

compare_utilities <- function(data_stan, sim_pop, est_utilities, nholdouts){
  if (sum(dim(sim_pop) != dim(est_utilities)) > 0){
    stop("Fatal Error: simulated and estimated utilities must be same size (rows and cols)")
  }
  plot(as.vector(sim_pop), as.vector(est_utilities), main = "Compare Utilities", xlab = "Simulated Population Utilities", ylab = "Estimated Utilities")
  plot(as.vector(cor(sim_pop)), as.vector(cor(est_utilities)), main = "Compare Correlation of Utilities", xlab = "Simulated Population Cor", ylab = "Estimated Correlation")
  plot(as.vector(cov(sim_pop)), as.vector(cov(est_utilities)), main = "Compare Covariance of Utilities", xlab = "Simulated Population Cov", ylab = "Estimated Covaraince")
  utility_stats <- list(
    MAE_utils = mean(abs(est_utilities - sim_pop)),
    cor_utils = cor(as.vector(est_utilities), as.vector(sim_pop)),
    mean_vars = data.frame(sim_pop = colMeans(sim_pop), est_utils = colMeans(est_utilities)),
    sd_vars = data.frame(sim_pop = apply(sim_pop, 2, sd), est_utils = apply(est_utilities, 2, sd)),
    sd_rows = data.frame(sim_pop = apply(sim_pop, 1, sd), est_utils = apply(est_utilities, 1, sd)),
    row_cor = data.frame(row_cor = sapply(1:nrow(sim_pop), function(i) cor(sim_pop[i,], est_utiities[i,])))
  )
  plot(utility_stats$mean_vars, main = "Compare Means of Utility Variables")
  plot(utility_stats$sd_vars, main = "Compare Std Deviations of Utility Variables")
  plot(utility_stats$sd_rows, main = "Compare Std Deviations of Each Respondent's Utilities")
  holdout_stats <- NULL; holdout_detail <- NULL
  if (nholdouts > 0){
    result <- do.call(rbind, lapply(1:nrow(sim_pop), function(i){
      MNL_holdout_compare_one(data_stan, sim_pop[i,], est_utilities[i,],
                              id = i,
                              hold_tasks = sample((1:data_stan$T)[data_stan$task_individual != i], nholdouts))
    }))
    holdout_stats <- list(
      hit_rate = sum(result$sim_choice * result$est_choice)/sum(result$sim_choice),
      rlh = exp(sum(result$sim_choice * log(result$est_pred))/sum(result$sim_choice)),
      MAE_probs = mean(abs(result$est_pred - result$sim_pred)),
      cor_probs = cor(result$sim_pred,result$est_pred)
    )
  }
  return(list(holdout_stats = holdout_stats, utility_stats = utility_stats, holdout_detail = result))
}

holdout_fit <- compare_utilities(data_stan, sim_pop, utilities, 30)
