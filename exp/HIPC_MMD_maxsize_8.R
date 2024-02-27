rm(list=ls())
set.seed(0)
setwd("/workdir2/solenne.gaucher/FlowCytometry")
#setwd("~/Desktop/Recherche/postdoc/FlowCyto/SimulationsV2")

# Load packages
library(doParallel)
library(foreach)
library(lubridate)

################################ Parameters ############################
cell_choice <- c(8)
theta_choice <-  c(0, exp(seq(log(0.001), log(0.3), length.out = 40)))
B_iterations <- 500
alpha_choice <- c(0, 0.003125, 0.00625, 0.0125, 0.025, 0.05, 0.06, 0.07)
max_size_test <- 10000
n_synth <- 10

################################ Load and pre-process HIPC data ############################
lab_names <- c("W2","FTV","IU","D54","O0","pM","pw")
data_df <- data.frame(matrix(ncol = 5+7, nrow = 0))
colnames(data_df) <- c("patient", "unique_patient", "replicate", "lab", "cell_type", 
                       "CCR7","CD4","CD45RA","CD3","HLADR","CD38","CD8")
unique_pat <- 1
for (lab in lab_names){
  for (pat in 1:3){
    for (rep in 1:3){
      if (file.exists(paste("data/HIPC/", lab,"_", 3*(pat-1) + rep, "_values.csv", sep= ""))){
        data_pat <- read.csv(paste("data/HIPC/", lab,"_", 3*(pat-1) + rep, "_values.csv", sep= ""))[,2:8]
        data_pat$cell_type <- read.csv(paste("data/HIPC/", lab,"_", 3*(pat-1) + rep, "_clust.csv", sep= ""))[, 2]
        data_pat$cell_type <- as.factor(data_pat$cell_type)
        data_pat$patient <- pat
        data_pat$patient <- as.factor(data_pat$patient)
        data_pat$unique_patient <- unique_pat
        data_pat$replicate <- rep
        data_pat$replicate <- as.factor(data_pat$replicate)
        data_pat$lab <- lab
        data_pat$lab <- as.factor(data_pat$lab)
        data_df <- rbind(data_df, data_pat)
        print(paste("lab : ", lab, ", pat : ", 3*(pat-1) + rep, ", sample size : ", nrow(data_pat)))
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
n.cores <- as.integer(parallel::detectCores()/6) - 1
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK",
  outfile = ""
)
doParallel::registerDoParallel(cl = my.cluster)

############################################################
print(paste("We begin the experiment, time is ", now()))
for (chosen_cell in cell_choice){
  data_chosen_cell <- data_df
  data_chosen_cell$cell_type <- 1*(data_chosen_cell$cell_type == chosen_cell)
  names(data_chosen_cell)[names(data_chosen_cell) == 'cell_type'] <- "Y"
  
  res_cell <- data.frame(matrix(nrow = 0, ncol = 5))
  names(res_cell) <- c("id", "theta",  "cell", "alpha", "power")
  
  res_cell <- foreach (id = levels(data_chosen_cell$unique_patient), .combine=rbind) %dopar%{
    source("utils/aux.R")
    source("utils/test_MMD.R")
    library(foreach)
    start <- lubridate::now()
    
    data_test <- data_chosen_cell[data_chosen_cell$unique_patient == id,]
    n_test <- min(max_size_test,nrow(data_test))
    data_train <- data_chosen_cell[data_chosen_cell$unique_patient != id & data_chosen_cell$Y == 0,]
    
    # Calibrate the test
    # Split data : half for computing bootstrap thresholds, half for validation
    id_thr <- sample(unique(data_train$unique_patient), floor(length(unique(data_train$unique_patient))/2))
    data_thr <- data_train[data_train$unique_patient %in% id_thr,]
    Xthr <- as.matrix(data_thr[,c("CCR7","CD4","CD45RA","CD3","HLADR","CD38","CD8")])
    data_val <- data_train[(!data_train$unique_patient %in% id_thr),]
    
    # Compute the quantiles of the test statistic for the different values of alpha 
    bootstrap_slice <- as.numeric(as.matrix(data_thr$unique_patient))
    qmmd <- quantile_MMDl_bootstrap(Xthr, n_test, B_iterations, alpha_choice, bootstrap_slice = bootstrap_slice, parallel = F)
    threshold <- qmmd$threshold
    sigma <- qmmd$sigma
    
    # Select a quantile best_alpha so that the test has empirical level 0 on validation set
    
    power <- rep(0, length(alpha_choice))
    for (id_val in unique(data_val$unique_patient)){
      Xval <- as.matrix(data_val[data_val$unique_patient == id_val,c("CCR7","CD4","CD45RA","CD3","HLADR","CD38","CD8")])
      Xval <- Xval[sample(nrow(Xval), n_test, replace = (nrow(Xval) < n_test)),]
      result  <- test_MMDl_bootstrap(Xthr, Xval, threshold, sigma)
      power <- power + result
    }
    if (power[1] > 0){
      best_alpha <- 0
    }else{
      best_alpha <- alpha_choice[max(which(power == 0))]
    }
    
    # Do the test!
    # Re-estimate the statistics and threshold
    Xtrain <- as.matrix(data_train[,c("CCR7","CD4","CD45RA","CD3","HLADR","CD38","CD8")])
    bootstrap_slice <- as.numeric(as.matrix(data_train$unique_patient))
    qmmd <- quantile_MMDl_bootstrap(Xtrain, n_test, B_iterations, alpha = best_alpha, bootstrap_slice = bootstrap_slice, parallel = F)
    threshold <- qmmd$threshold
    sigma <- qmmd$sigma
    
    ######################################
    # To create synthetic test sets for different prevalence theta, split cells 0 and 1
    Xtest_0 <- as.matrix(data_test[data_test$Y == 0,c("CCR7","CD4","CD45RA","CD3","HLADR","CD38","CD8")])
    Xtest_1 <- as.matrix(data_test[data_test$Y == 1,c("CCR7","CD4","CD45RA","CD3","HLADR","CD38","CD8")])
    res_id <- foreach (theta = theta_choice, .combine=rbind) %do% {
      reject <- 0
      for (r in 1:n_synth){
        if (as.integer(theta*n_test) > 0){
          Xtest <- rbind(Xtest_0[sample(nrow(Xtest_0), n_test - as.integer(theta*n_test), 
                                        replace = (nrow(Xtest_0) < n_test - as.integer(theta*n_test))), ], 
                         Xtest_1[sample(nrow(Xtest_1), as.integer(theta*n_test),
                                        replace = (nrow(Xtest_1) < as.integer(theta*n_test))), ])
        }else{
          Xtest <- Xtest_0[sample(nrow(Xtest_0), n_test, replace = TRUE), ]
        }
        reject <- reject + test_MMDl_bootstrap(Xtrain, Xtest, threshold, sigma)
      }
      return(data.frame(id = id, theta = theta, cell = chosen_cell, best_alpha = best_alpha, power = reject/n_synth, max_size_test = max_size_test))
    }
    saveRDS(res_id, paste("results/HIPC_MMD_maxsize/", chosen_cell, "_", id, ".rds", sep = ""))
    print(paste("Computation done for cell :", chosen_cell, ", with id ", id))
    print(lubridate::now() - start)
    print("")
    print("")
    return(res_id)
  }
  print(paste("results/HIPC_MMD_maxsize/", chosen_cell, "tot.rds", sep = ""))
  saveRDS(res_cell, paste("results/HIPC_MMD_maxsize/", chosen_cell, "tot.rds", sep = ""))
}

parallel::stopCluster(cl = my.cluster)