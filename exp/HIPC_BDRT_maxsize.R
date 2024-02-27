rm(list=ls())
set.seed(0)
setwd("/workdir2/solenne.gaucher/FlowCytometry")
#setwd("~/Desktop/Recherche/postdoc/FlowCyto/SimulationsV2")

# Load packages
library(doParallel)
library(partykit)
library(foreach)
library(lubridate)

################################ Parameters ############################
n_sample_part <- 12 # 12 samples are used for choosing the partition
cell_choice <- c(2, 6, 8, 9)
theta_choice <-  c(0, exp(seq(log(0.001), log(0.3), length.out = 40)))
n_synth <- 10 # number synthetic datasets created from each sample at test time
B_iterations <- 500
alpha_choice <- c(0, 0.003125, 0.00625, 0.0125, 0.025, 0.05)
max_size_test <- 10000
alpha <- 0.5

################################ Load and pre-process HIPC data ############################
lab_names <- c("W2","FTV","IU","D54","O0","pM","pw")
data_df <- data.frame(matrix(ncol = 5+7, nrow = 0))
colnames(data_df) <- c("patient", "unique_patient", "replicate", "lab", "cell_type", 
                       "CCR7","CD4","CD45RA","CD3","HLADR","CD38","CD8")
unique_pat <- 1
for (lab in lab_names){
  for (pat in 1:3){
    for (replicate in 1:3){
      if (file.exists(paste("data/HIPC/", lab,"_", 3*(pat-1) + replicate, "_values.csv", sep= ""))){
        data_pat <- read.csv(paste("data/HIPC/", lab,"_", 3*(pat-1) + replicate, "_values.csv", sep= ""))[,2:8]
        data_pat$cell_type <- read.csv(paste("data/HIPC/", lab,"_", 3*(pat-1) + replicate, "_clust.csv", sep= ""))[, 2]
        data_pat$cell_type <- as.factor(data_pat$cell_type)
        data_pat$patient <- pat
        data_pat$patient <- as.factor(data_pat$patient)
        data_pat$unique_patient <- unique_pat
        data_pat$replicate <- replicate
        data_pat$replicate <- as.factor(data_pat$replicate)
        data_pat$lab <- lab
        data_pat$lab <- as.factor(data_pat$lab)
        data_df <- rbind(data_df, data_pat)
        print(paste("lab : ", lab, ", pat : ", 3*(pat-1) + replicate, ", sample size : ", nrow(data_pat)))
        print(paste("data/HIPC/", lab,"_", 3*(pat-1) + replicate, "_values.csv", sep= ""))
        unique_pat <- unique_pat + 1
      }
    }
  }
}

### Scale data to [0,1]**7 ###
data_df$unique_patient <- as.factor(data_df$unique_patient)
data_df[, c("CCR7","CD4","CD45RA","CD3","HLADR","CD38","CD8")] <- apply(data_df[,c("CCR7","CD4","CD45RA","CD3","HLADR","CD38","CD8")], MARGIN = 2,
                                                                        FUN = function(X) (X - mean(X)))
data_df[, c("CCR7","CD4","CD45RA","CD3","HLADR","CD38","CD8")] <- asinh(data_df[, c("CCR7","CD4","CD45RA","CD3","HLADR","CD38","CD8")])
data_df[, c("CCR7","CD4","CD45RA","CD3","HLADR","CD38","CD8")] <- apply(data_df[,c("CCR7","CD4","CD45RA","CD3","HLADR","CD38","CD8")], MARGIN = 2, 
                                                                        FUN = function(X) (X - min(X))/(range(X)[2]- range(X)[1]))

################################ For parallelization ############################
n.cores <- as.integer(parallel::detectCores()/3) - 1
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK",
  outfile = ""
)
doParallel::registerDoParallel(cl = my.cluster)

############################################################
print(paste("We begin the experiment, time is ", now()))
for (chosen_cell in cell_choice){
  # label Y = 1 if cell == chosen_cell, 0 otherwise
  data_chosen_cell <- data_df
  data_chosen_cell$cell_type <- 1*(data_chosen_cell$cell_type == chosen_cell)
  names(data_chosen_cell)[names(data_chosen_cell) == 'cell_type'] <- "Y"
  
  res_cell <- foreach (id = unique(data_chosen_cell$unique_patient), .combine=rbind) %dopar%{
    # conduct the experiment one test sample at a time
    source("utils/user_defined_split.R")
    source("utils/aux.R")
    source("utils/test_dr.R")
    source("utils/part_choice.R")
    library(foreach)
    
    start <- lubridate::now()
    
    data_test <- data_chosen_cell[data_chosen_cell$unique_patient == id,]
    n_test <- min(max_size_test,nrow(data_test))
    data_train <- data_chosen_cell[data_chosen_cell$unique_patient != id,]
    
    # Split data_train for partition and estimation
    id_part <- sample(unique(data_train$unique_patient), n_sample_part)
    data_part <- data_train[data_train$unique_patient %in% id_part,]
    data_est <- data_train[!data_train$unique_patient %in% id_part, ]
    
    n1_part <- sum(data_part$Y)
    n0_part <- nrow(data_part) - n1_part
    n1_est <- sum(data_est$Y)
    n0_est <- nrow(data_est) - n1_est
    
    # Define the parameters for growing the partition sequence
    Kmax <- 15
    u <- log(4*Kmax/alpha)
    t <- log(2/alpha)
    epsilon0 = max(3*u/n0_part, t/n_test)
    epsilon1 = max(sqrt(3*u/n1_part), epsilon0)
    parms_part = list(epsilon0 = epsilon0, epsilon1 = epsilon1, n0 = n0_part, n1 = n1_part)
    
    # Choose the partition sequence using respectively the density ratio criteria and the gini id
    full_tree <- rpart::rpart(Y ~ CCR7 + CD4 + CD45RA + CD3 + HLADR + CD38 + CD8, data = data_part,
                              method = method_dr(), parms = parms_part, cp = 0, maxdepth = 10, minbucket = sqrt(n1_part))
    pruned_tree <- prune_dr(full_tree, data_part, parms_part, Kmin = 10)
    partition <- as.party(pruned_tree$pruned_tree)
    
    # Calibrate the test
    # Split data: half for computing bootstrap thresholds, half for validation
    id_thr <- sample(unique(data_est$unique_patient), floor(length(unique(data_est$unique_patient))/2))
    data_thr <- data_est[data_est$unique_patient %in% id_thr,]
    data_val <- data_est[!data_est$unique_patient %in% id_thr, ]
    
    # Compute the quantiles of the test statistic for the different values ofalpha_choice_choice
    u <- log(4*Kmax/alpha)
    t <- log(2/alpha)
    epsilon0 = max(3*u/n0_est, t/n_test)
    epsilon1 = max(sqrt(3*u/n1_est),epsilon0)
    parms_est = list(epsilon0 = epsilon0, epsilon1 = epsilon1, n0 = n0_est, n1 = n1_est,
                     u = u, t = t, n_test = n_test)
    threshold <- quantile_dr_bootstrap(data_thr, partition, parms_est, B_iterations, alpha_choice,
                                       bootstrap_slice = data_thr$unique_patient)
    
    # Select a quantile best_alpha so that the test has empirical level 0 on validation set
    power <- rep(0, length(alpha_choice))
    for (id_val in unique(data_val$unique_patient)){
      data_one_val <- data_val[data_val$unique_patient == id_val & data_val$Y == 0,]
      data_one_val <- data_one_val[sample(nrow(data_one_val), n_test, replace = (nrow(data_one_val) < n_test)),]
      power <- power + test_dr_bootstrap(data_thr, data_one_val, partition, parms_est, threshold)
    }
    if (power[1] > 0){
      best_alpha <- 0
    }else{
      best_alpha <- alpha_choice[max(which(power == 0))]
    }
    
    # Do the test!
    # Re-estimate the statistics and threshold
    threshold <- quantile_dr_bootstrap(data_est, partition, parms_est, B_iterations, best_alpha,
                                       bootstrap_slice = data_est$unique_patient)
    data_test_0 <- data_test[data_test$Y == 0, ]
    data_test_1 <- data_test[data_test$Y == 1, ]
    
    res_id <- foreach (theta = theta_choice, .combine=rbind) %do% {
      reject <- 0
      for (r in 1:n_synth){
        if (as.integer(theta*n_test) > 0){
          data_test <- rbind(data_test_0[sample(nrow(data_test_0), n_test - as.integer(theta*n_test), 
                                                replace = (nrow(data_test_0) < n_test - as.integer(theta*n_test))), ], 
                             data_test_1[sample(nrow(data_test_1), as.integer(theta*n_test),
                                                replace = (nrow(data_test_1) < as.integer(theta*n_test))), ])
        }else{
          data_test <- data_test_0[sample(nrow(data_test_0), n_test, replace = TRUE), ]
        }
        reject <- reject + test_dr_bootstrap(data_est, data_test, partition, parms_est, threshold)
      }
      return(data.frame(id = id, theta = theta, cell = chosen_cell, best_alpha = best_alpha, power = reject/n_synth, max_size_test = max_size_test))
    }
    
    saveRDS(res_id, paste("results/HIPC_BDRT_maxsize/", chosen_cell, "_", id, ".rds", sep = ""))
    print(paste("Computation done for cell :", chosen_cell, ", with id ", id))
    print(lubridate::now() - start)
    print("")
    print("")
    return(res_id)
  }
  
  saveRDS(res_cell, paste("results/HIPC_BDRT_maxsize/", chosen_cell, "tot.rds", sep = ""))
}
parallel::stopCluster(cl = my.cluster)

