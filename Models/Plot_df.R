draws <- x
library(posterior)
saveRDS(draws_beta, "draws_beta.rds")
draws_test <- readRDS("draws_beta.rds")

plot_draws_df <- function(draws, vnames = NULL, ylab = "Draw",
                          chain_colors = rep(c("red","blue","green","black"),2)){
  # Assume draws is a list of data frames
  # Each data frame is draws for one chain of results with rows for iterations 
  # We plot each column
  if (is.null(vnames)) vnames <- colnames(draws[[1]]) # take names from first data frame
  fit_stats <- data.frame(
    variable = vnames,
    mean = NA,
    sd =  NA,
    rhat = NA,
    ESS = NA
  )
  for (i in 1:ncol(draws[[1]])){
    x <- sapply(1:length(draws), function(chain){
      draws[[chain]][,i]     
    })
    fit_stats$mean[i] <- round(mean(x), 2)
    fit_stats$sd[i] <- round(sd(x),2)
    fit_stats$rhat[i] <- round(rhat(x),2)
    fit_stats$ESS[i] <- round(ess_basic(x),1)
    plot(x[,1], type = "l", col = chain_colors[1], ylim = c(min(x), max(x)),
         xlab = "Sample Iteration", ylab = "Draw",
         main = paste(vnames[i],
                      "| rhat = ", round(rhat(x),2),
                      "| ESS = ", round(ess_basic(x),1)
         ))
    if (ncol(x) > 1){
      for (chain in 2:ncol(x)){
        lines(x[,chain], type = "l", col = chain_colors[chain])
      }
    }
  } # end for
  return(fit_stats)
}

fit_kl <- plot_draws_df(draws_test$post_warmup_draws)
