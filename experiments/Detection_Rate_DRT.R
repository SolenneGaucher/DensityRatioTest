rm(list=ls())
set.seed(0)
setwd("/workdir2/solenne.gaucher/FlowCytometry")
#setwd("~/Desktop/Recherche/postdoc/FlowCyto/SimulationsV2")

# Load packages
library(doParallel)
library(partykit)
library(foreach)
library(lubridate)

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
# n_train = n0_train + n1_train
n_train_choice <- as.integer(exp(seq(log(1000), log(1000000), length.out = 10)))+1 # values for n_train
p0 <- 0.7
theta_choice <-  c(0, exp(seq(log(0.0003), log(0.3), length.out = 30))) # values for theta
n_sim <- 100 # number of simulations
B_iterations <- 300 # number of repetitions for the bootstrap quantile
K_max <- 15 # maximal number of bins
alpha <- 0.05

################################ For parallelization ############################
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK",
  outfile = ""
)
doParallel::registerDoParallel(cl = my.cluster)

########################### Simulations ########################### 
print(paste("We begin the experiment, time is ", now()))
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
      
      n_test <- as.integer(n_train/10)
      
      # Define the parameters for growing the partition sequence
      Kmax <- 15
      u <- log(4*Kmax/alpha)
      t <- log(2/alpha)
      epsilon0 = max(3*u/n0_part, t/n_test)
      epsilon1 = max(sqrt(3*u/n1_part),epsilon0)
      parms_part = list(epsilon0 = epsilon0, epsilon1 = epsilon1, n0 = n0_part, n1 = n1_part)
      
      # Choose the partition sequence using respectively the density ratio criteria and the gini index
      full_tree <- rpart::rpart(Y ~ X1 + X2, data = data_part,
                                method = method_dr(), parms = parms_part, cp = 0, maxdepth = 10, minbucket = sqrt(n1_part))
      
      # Define the parameters for choosing the best partition
      u <- log(4*Kmax/alpha)
      t <- log(2/alpha)
      epsilon0 = max(3*u/n0_est, t/n_test)
      epsilon1 = max(sqrt(3*u/n1_est),epsilon0)
      parms_est = list(epsilon0 = epsilon0, epsilon1 = epsilon1, n0 = n0_est, n1 = n1_est)
      
      # Choose the best partition in the sequence
      pruned_tree <- prune_dr(full_tree, data_est, parms_est, Kmin = 3)
      partition <- as.party(pruned_tree$pruned_tree)
      
      # Set the parameters for the test
      K <- pruned_tree$K
      u <- log(4*K/alpha)
      t <- log(2/alpha)
      epsilon0 = max(3*u/n0_est, t/n_test)
      epsilon1 = max(sqrt(3*u/n1_est),epsilon0)
      parms_est = list(epsilon0 = epsilon0, epsilon1 = epsilon1, 
                       n0 = n0_est, n1 = n1_est, n_test = n_test,
                       u = u, t = t, K = K)
      
      # Conduct the test
      res <- as.data.frame(matrix(nrow = 0, ncol = 5))
      for (theta in theta_choice){
        n1_test <- as.integer(theta*n_test)
        n0_test <- n_test - n1_test
        if (n1_test >= 1){
          Xtest <- rbind(rf0(n0_test, setting),rf1(n1_test, setting))
        }else{
          Xtest <- rf0(n0_test, setting)
        }
        data_test <- data.frame(Xtest)
        res <- rbind(res, c("DRT", setting, n_train, theta, 
                            test_dr(data_est, data_test, parms_est, partition)))
      }
      names(res) <- c("method" , "setting", "n_train", "theta", "result")
      return(res)
    }
    print(paste("setting = ", setting,", n_train = ", n_train, ", time is ", now()))
    saveRDS(res_exp, paste("results/Detection_rate_DRT/setting_",setting, "_ntrain", n_train,".rds", sep = ""))
  }
}
parallel::stopCluster(cl = my.cluster)
