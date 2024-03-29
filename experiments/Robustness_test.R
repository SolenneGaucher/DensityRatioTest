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

################################## Distributions ##############################
s <- matrix(c(1/100, 0, 0, 1/100), ncol = 2, nrow = 2)
# Define distribution f0

m0a <- c(3/10, 6/10)
m0b <- c(6/10, 3/10)
rf0 <- function(n, p){
  return(rbind(tmvtnorm::rtmvnorm(as.integer(n*p), m0a, s, lower = c(0,0), upper = c(1,1)),
               tmvtnorm::rtmvnorm(n - as.integer(n*p), m0b, s, lower = c(0,0), upper = c(1,1))))
}

# Distribution f1
m1 <- c(6/10, 6/10)
rf1 <- function(n){
  return(tmvtnorm::rtmvnorm(n, m1, s, lower = c(0,0), upper = c(1,1)))}

################################ Choose parameters for this simulation ####################################################
n_sim <- 100 # number of simulation
n_train <- 100000
p0 <- 0.7 # n0/ntrain 
p_choice <- c(0.5, 0.6, 0.7, 0.8, 0.9) 
exp_choice <- c("I", "II")
# in exp I, we fix p = 0.5 in training data and vary it in test date
# in exp II we do the opposite
theta_choice <- c(0, 0.015)


################################ For parallelization ############################
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

########################### Simulations ########################### 
print(paste("We begin the experiment, time is ", now()))
for (expe in exp_choice){
  for (p in p_choice){
    res_exp <- foreach(i = 1:n_sim, .combine=rbind) %dopar% {
      source("utils/user_defined_split.R")
      source("utils/aux.R")
      source("utils/test_dr.R")
      source("utils/part_choice.R")
      library(lubridate)
      
      # Simulate data for X_part^0, X_part^1, X_est^0, X_est^1
      n_part <- as.integer(n_train/2)
      n0_part <- as.integer(p0*n_train/2)
      n1_part <- n_part - n0_part
      
      n_est <- as.integer(n_train/2)
      n0_est <- as.integer(p0*n_train/2)
      n1_est <- n_est - n0_est
      
      print(paste("sim ", i,", flag 0"))
      if (expe == "II"){
        Xpart <- rbind(rf0(n0_part, p),rf1(n1_part))
        Xest <- rbind(rf0(n0_est, p),rf1(n1_est))
      }else{
        Xpart <- rbind(rf0(n0_part, 0.5),rf1(n1_part))
        Xest <- rbind(rf0(n0_est, 0.5),rf1(n1_est))
      }
      
      data_part <- data.frame(Xpart)
      data_part$Y <- c(rep(0, n0_part), rep(1, n1_part))
      data_est <- data.frame(Xest)
      data_est$Y <- c(rep(0, n0_est), rep(1, n1_est))
      
      # Define the parameters 
      Kmax <- 15
      u_part <- log(4*Kmax/0.05)
      parms_part = list(epsilon0 = 3*u_part/n0_part, epsilon1 = max(sqrt(3*u_part/n1_part),3*u_part/n0_part), n0 = n0_part, n1 = n1_part)
      
      u_est <- log(4*Kmax/0.05)
      parms_est = list(epsilon0 = 3*u_est/n0_est, epsilon1 = max(sqrt(3*u_est/n1_est),9*u_est/n0_est),
                       n0 = n0_est, n1 = n1_est, u = u_est, t = log(2/0.05))
      
      # Choose the partition sequence using respectively the density ratio criteria and the gini index
      # full_tree <- rpart::rpart(Y ~ X1 + X2, data = data_part,
      #                             method = method_dr(), parms = parms_part, cp = 0, maxdepth = 10, minbucket = sqrt(n1_part))
      
      full_tree <- tryCatch({rpart::rpart(Y ~ X1 + X2, data = data_part, method = method_dr(), parms = parms_part, cp = 0, maxdepth = 10, minbucket = sqrt(n1_part))},
               error = function(e) {print(paste("error at full_tree, sim is ", i))})
      # Choose the best partition in the sequence
      pruned_tree <- tryCatch({prune_dr(full_tree, data_est, parms_est, verbose = F, Kmin = 3)},
                            error = function(e) {print(paste("error at pruned_tree, sim is ", i))})
      # pruned_tree <- prune_dr(full_tree, data_est, parms_est, verbose = F, Kmin = 3)
      # partition <- as.party(pruned_tree$pruned_tree)
      partition <- tryCatch({as.party(pruned_tree$pruned_tree)},
                            error = function(e) {print(paste("error at party, sim is ", i))})
      
      # Estimate the signal for each partition
      K <- pruned_tree$K
      u_est <- log(4*K/0.05)
      parms_est = list(epsilon0 = 3*u_est/n0_est, epsilon1 = max(sqrt(3*u_est/n1_est),9*u_est/n0_est),
                       n0 = n0_est, n1 = n1_est, u = u_est, t = log(2/0.05), K = K)
      stats <- compute_h_partition(partition, data_est, parms_est)
      
      # Conduct the test
      n_test <- as.integer(n_train/10)
      res <- as.data.frame(matrix(nrow = 0, ncol = 5))
      for (theta in theta_choice){
        n1_test <- as.integer(theta*n_test)
        n0_test <- n_test - n1_test
        if (n1_test >= 1){
          if (expe == "I"){
            Xtest <- rbind(rf0(n0_part, p),rf1(n1_part))
          }else{
            Xtest <- rbind(rf0(n0_part, 0.5),rf1(n1_part))
          }
        }else{
          if (expe == "I"){
            Xtest <- rf0(n0_part, p)
          }else{
            Xtest <- rf0(n0_part, 0.5)
          }
        }
        data_test <- data.frame(Xtest)
        test_res <- tryCatch({test_dr(data_est, data_test, parms_est, partition)},
                                                            error = function(e) {print(paste("error at test, sim is ", i, ", theta is ", theta))})
        res <- rbind(res, c("DRT", expe, theta, test_res, p))
        # res <- rbind(res, c("DRT", expe, theta, test_dr(data_est, data_test, parms_est, partition)))
      }
      names(res) <- c("method", "expe", "theta", "result", "p")
      return(res)
    }
    print(paste("exp = ", expe,", n_train = ", n_train, ", p = ", p, ", time is ", now()))
    saveRDS(res_exp, paste("results/Robustness/expe_",expe,"_p_", p, ".rds", sep = ""))
  }
}
parallel::stopCluster(cl = my.cluster)

