cov <- data_stan$prior_cov
I <- 300
draw_alpha <- FALSE
draw_cov <- TRUE
con_sign <- data_stan$con_sign
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