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
n_train_choice <- as.integer(exp(seq(log(1000), log(1000000), length.out = 10)))[1:7]+1 # values for n_train
p0 <- 0.7
theta_choice <-  c(0, exp(seq(log(0.0003), log(0.3), length.out = 30))) # values for theta
n_sim <- 100 # number of simulations
B_iterations <- 300 # number of repetitions for the bootstrap quantile
alpha <- 0.5

################################ For parallelization ############################
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

########################### Simulations ########################### 
print(paste("We begin the experiment, time is ", now()))
for (setting in c("A", "B", "C")){
  for (n_train in n_train_choice){
    res_exp <- foreach(i = 1:n_sim, .combine=rbind) %dopar% {
      source("utils/user_defined_split.R")
      source("utils/aux.R")
      source("utils/test_MMD.R")
      source("utils/part_choice.R")
      
      # Simulate data for X^0
      n0_train <- as.integer(p0*n_train)
      Xtrain <- rf0(n0_train, setting)
      
      # Compute the threshold for the bootstrap MMD
      n_test <- as.integer(n_train/10)
      stats_mmd <- quantile_MMDl_bootstrap(Xtrain, n_test, B_iterations)
      threshold <- stats_mmd$threshold
      sigma <- stats_mmd$sigma
      
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
        Xtest <- as.matrix(Xtest)
        res <- rbind(res, c("MMD", setting, n_train, theta, test_MMDl_bootstrap(Xtrain, Xtest, threshold, sigma)))
      }
      names(res) <- c("method", "setting", "n_train", "theta", "result")
      return(res)
    }
    print(paste("setting = ", setting,", n_train = ", n_train, ", time is ", now()))
    saveRDS(res_exp, paste("results/Detection_rate_MMD/setting_",setting, "_ntrain", n_train,".rds", sep = ""))
  }
}
parallel::stopCluster(cl = my.cluster)
