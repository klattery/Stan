library(sqldf)

dir <- "C:/Users/K.Lattery/SKIM/LinkedIn - 4. Analysis/3. HB_Stan/Kevin"
data_in <- read.csv(file.path(dir,"Data_Kai.csv"), check.names = FALSE)
prod_recode <- read.csv(file.path(dir,"Product_Recodev3.csv"))
prod_vars <- prod_recode[data_in$Product,]
data_merge <- cbind(data_in, prod_vars)
data_merge$price1 <- data_merge$Price * (data_merge$Product == 1)
data_merge$price2 <- data_merge$Price * (data_merge$Product == 2)
data_merge$price3 <- data_merge$Price * (data_merge$Product == 3)
data_merge$price4 <- data_merge$Price * (data_merge$Product == 4)
data_merge$price5 <- data_merge$Price * (data_merge$Product == 5)
data_merge$price_posts <- data_merge$Price * (data_merge$Product == 6)

########################
# Licenses Products 1- 5
data_license <- data_merge

num_license_task <- sqldf("select sys_Respnum, Task, sum(Response) as Num_Licenses
                          from data_license where Product < 6 group by sys_Respnum, Task")
num_license_max <- sqldf("select sys_Respnum, max(Num_Licenses) as Max_Licenses
                         from num_license_task group by sys_Respnum")
num_license <- merge(num_license_task, num_license_max)
num_license$prob1 <- num_license$Num_Licenses/num_license$Max_Licenses
num_license$prob1[is.na(num_license$prob1)] <- 0
num_license$prob2 <- 1 - num_license$prob1
num_license$multresp <- num_license$prob1/num_license$Num_Licenses
num_license$multresp[is.na(num_license$multresp)] <- 0

# Add num_license data
data_license <- merge(data_license, num_license)
data_license$dep <- data_license$Response * data_license$multresp
data_license$dep[data_license$Product == 6] <- 0

# Add None rows
none <- data.frame("sys_RespNum" = num_license$sys_RespNum,
                   "Task" = num_license$Task,
                   "Concept" = 99,
                   dep = num_license$prob2)
data_license_wnone <- rbind(
  data_license,
  data.frame(c(none, sapply(setdiff(names(data_license), names(none)), function(x) 0)), check.names = FALSE)
)
data_license_wnone <- data_license_wnone[order(data_license_wnone$sys_RespNum,
                                               data_license_wnone$Task,
                                               data_license_wnone$Concept),]
data_license_wnone$None1_5 <- as.numeric(data_license_wnone$Product == 0)


# Add in assumptions of Premium
prem <- (data_license_wnone$Product == 5)
data_license_wnone$`Network Access`[prem] <- 2 # Unlimited
data_license_wnone$`Spotlights/Open to Work Filter`[prem] <- 1 # No
data_license_wnone$`InMail Messaging`[prem] <- 1.5 # 15 is halfway between None (1) and 30 (2)
data_license_wnone$`Bulk InMails/Messaging`[prem] <- 1 # No
data_license_wnone$`Job posts`[prem] <- 1 # Only 1 post

# Add level to job_posts attribute based on whether ala carte option exists
data_license_wnone$Job_Posts_r <- data_license_wnone$`Job posts`
ksplit <- split(data_license_wnone, list(data_license_wnone$Task, data_license_wnone$sys_RespNum), drop = TRUE)
kdata_list <- lapply(ksplit, function(x){
  result <- x
  if (!(6 %in% result$Product)){
    result$Job_Posts_r[result$Job_Posts_r == 2] <- 3 # Recode level 2 to level 3 
  } 
  return(result)
})
data_license_final <- do.call(rbind, kdata_list)

# Move dep to end and remove a few unneeded vars
dep <- data_license_final$dep
data_license_final <- data_license_final[order(data_license_final$sys_RespNum, data_license_final$Task, data_license_final$Concept),
                                         !(colnames(data_license_final) %in% c("prob1","prob2","multresp","Job posts", "dep"))]
data_license_final$task_type <- 1
data_license_final$dep <- dep

######################################
#  Get Optional Job Posts Max
num_posts_max <- sqldf("select sys_Respnum, max(Response) as Max_Posts
                         from data_merge where data_merge.Product = 6 group by sys_Respnum")

write.csv(data_license_final, file.path(dir, "Comb_model","data_license_final_v4.csv"), row.names = FALSE)
write.csv(num_license_max, file.path(dir, "Comb_model", "num_license_max.csv"), row.names = FALSE)
write.csv(num_posts_max, file.path(dir,"Comb_model","num_posts_max.csv"), row.names = FALSE)


