cat("Number of concepts per task")
as.data.frame(table(data_stan$end - data_stan$start + 1, dnn = "Concepts in Task"), responseName = "Num Tasks")
