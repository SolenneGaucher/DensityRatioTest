rm(list=ls())
set.seed(0)
setwd("/workdir2/solenne.gaucher/FlowCytometry")
#setwd("~/Desktop/Recherche/postdoc/FlowCyto/SimulationsV2")

# Load packages
library(doParallel)
library(partykit)
library(foreach)
library(doParallel)
library(lubridate)

################################ For parallelization ############################
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

################################## Distributions ##############################
s <- matrix(c(1/100, 0, 0, 1/100), ncol = 2, nrow = 2)
# Define distribution f0

rf0 <- function(n, setting){
  if (setting == "A"){
    m0 <- c(3/10, 3/10)
  }else if (setting == "B"){
    m0 <- c(4/10, 4/10)
  }else{
    m0 <- c(4/10, 4/10)
  }
  return(tmvtnorm::rtmvnorm(n, m0, s, lower = c(0,0), upper = c(1,1)))
}

# Define distribution f1
rf1 <- function(n, setting){
  if (setting == "A"){
    m0 <- c(7/10, 7/10)
  }else if (setting == "B"){
    m0 <- c(6/10, 6/10)
  }else{
    m0 <- c(5/10, 5/10)
  }
  return(tmvtnorm::rtmvnorm(n, m0, s, lower = c(0,0), upper = c(1,1)))
}

################################ Choose parameters for this simulation ####################################################
n_sim <- 500 # number of simulation
n_train_choice <- c(1000, 3000, 10000, 30000, 100000, 300000, 1000000) # size of training sample X^0 + size of training sample X^1
p0 <- 0.7 # n0/(n0+n1)


########################### Simulations ########################### 
for (setting in c("A", "B", "C")){
  for (n_train in n_train_choice){
    res_exp <- foreach(i = 1:n_sim, .combine=rbind) %dopar% {
      source("utils/user_defined_split.R")
      source("utils/aux.R")
      source("utils/test_dr.R")
      source("utils/part_choice.R")
      
      # Simulate data for X_part^0, X_part^1, X_est^0, X_est^1
      n0_part <- as.integer(p0*n_train/2)
      n1_part <- as.integer((1-p0)*n_train/2)
      n_part <- n0_part + n1_part
      Xpart <- rbind(rf0(n0_part, setting),rf1(n1_part, setting))
      data_part <- data.frame(Xpart)
      data_part$Y <- c(rep(0, n0_part), rep(1, n1_part))
      
      n0_est <- as.integer(p0*n_train/2)
      n1_est <- as.integer((1-p0)*n_train/2)
      n_est <- n0_est + n1_est
      Xest <- rbind(rf0(n0_est, setting),rf1(n1_est, setting))
      data_est <- data.frame(Xest)
      data_est$Y <- c(rep(0, n0_est), rep(1, n1_est))
      
      # Define the parameters 
      Kmax <- 15
      u_part <- log(4*Kmax/0.05)
      parms_part = list(epsilon0 = 3*u_part/n0_part, epsilon1 = max(sqrt(3*u_part/n1_part),3*u_part/n0_part), n0 = n0_part, n1 = n1_part)
      
      u_est <- log(4*Kmax/0.05)
      parms_est = list(epsilon0 = 3*u_est/n0_est, epsilon1 = max(sqrt(3*u_est/n1_est),9*u_est/n0_est), n0 = n0_est, n1 = n1_est)
      
      # Choose the partition sequence using respectively the density ratio criteria and the gini index
      full_tree_dr <- rpart::rpart(Y ~ X1 + X2, data = data_part,
                                   method = method_dr(), parms = parms_part, cp = 0, maxdepth = 10, minbucket = sqrt(n1_part))
      full_tree_gini <- rpart::rpart(Y ~ X1 + X2, data = data_part,
                                     method = "class", cp = 0, maxdepth = 10, minbucket = sqrt(n1_part))
      
      # Choose the best partition in the sequence
      pruned_tree_dr <- as.party(prune_dr(full_tree_dr, data_est, parms_est, verbose = F, Kmin = 3)$pruned_tree)
      pruned_tree_gini<- as.party(prune_dr(full_tree_gini, data_est, parms_est, verbose = F, Kmin = 3)$pruned_tree)
      
      # Estimate the signal for each partition
      stats_dr <- compute_h_partition(pruned_tree_dr, data_est, parms_est)
      sigma_dr <- sum((stats_dr$h1k/stats_dr$h0k-1)**2 * stats_dr$h0k)
      stats_gini <- compute_h_partition(pruned_tree_gini, data_est, parms_est)
      sigma_gini <- sum((stats_gini$h1k/stats_gini$h0k-1)**2 * stats_gini$h0k)
      
      # return the results
      res <- rbind(c(setting, n_train, "dr", "no", sigma_dr), 
                   c(setting, n_train, "gini", "no", sigma_gini))
      res <- as.data.frame(res)
      names(res) <- c("setting", "n_train", "criterion", "Kmin", "sigma")
      
      if (setting == "C"){
        # Additional simulation with Kmin = log(n_train)
        # Choose the best partition in the sequence
        pruned_tree_dr <- as.party(prune_dr(full_tree_dr, data_est, parms_est, verbose = F, Kmin = log(n_train))$pruned_tree)
        pruned_tree_gini<- as.party(prune_dr(full_tree_gini, data_est, parms_est, verbose = F, Kmin = log(n_train))$pruned_tree)
        
        # Estimate the signal for each partition
        stats_dr <- compute_h_partition(pruned_tree_dr, data_est, parms_est)
        sigma_dr <- sum((stats_dr$h1k/stats_dr$h0k-1)**2 * stats_dr$h0k)
        stats_gini <- compute_h_partition(pruned_tree_gini, data_est, parms_est)
        sigma_gini <- sum((stats_gini$h1k/stats_gini$h0k-1)**2 * stats_gini$h0k)
        
        # return the results
        res <- rbind(res,
                     c(setting, n_train, "dr", "log(n)", sigma_dr), 
                     c(setting, n_train, "gini", "log(n)", sigma_gini))
        names(res) <- c("setting", "n_train", "criterion", "Kmin", "sigma")
      }
      return(res)
    }
    print(paste("setting = ", setting,", n_train = ", n_train, ", time is ", now()))
    saveRDS(res_exp, paste("results/DROP_ABC/setting_",setting, "_ntrain", n_train,".rds", sep = ""))
  }
}
parallel::stopCluster(cl = my.cluster)
