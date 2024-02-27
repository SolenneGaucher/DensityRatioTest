rm(list=ls())
set.seed(0)
setwd("/workdir2/solenne.gaucher/FlowCytometry")
#setwd("~/Desktop/Recherche/postdoc/FlowCyto/SimulationsV2")

# Load packages
library(doParallel)
library(foreach)
library(lubridate)

################################ Parameters ############################
cell_choice <- c(2, 6, 8, 9)
theta_choice <-  c(0, exp(seq(log(0.001), log(0.3), length.out = 40)))
rep <- 10 # number synthetic datasets created from each sample at test time
B_bootstrap <- 500
alpha <- c(0, 0.003125, 0.00625, 0.0125, 0.025, 0.05, 0.06, 0.07)

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
        print(paste("data/HIPC/", lab,"_", 3*(pat-1) + rep, "_values.csv", sep= ""))
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
n.cores <- parallel::detectCores() - 1
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
    source("aux.R")
    source("test_MMD.R")
    library(lubridate)
    
    start_time <- Sys.time()
    data_test <- data_chosen_cell[data_chosen_cell$unique_patient == id,]
    data_train <- data_chosen_cell[data_chosen_cell$unique_patient != id & data_chosen_cell$Y == 0,]
    n_test <- nrow(data_test)
    
    # Select a quantile best_alpha so that the test has empirical level 0.05 on validation set
    # Split data : half for computing bootstrap thresholds, half for validation
    id_thr <- sample(unique(data_train$unique_patient), floor(length(unique(data_train$unique_patient))/2))
    data_thr <- data_train[data_train$unique_patient %in% id_thr,]
    Xthr <- as.matrix(data_thr[,c("CCR7","CD4","CD45RA","CD3","HLADR","CD38","CD8")])
    data_val <- data_train[(!data_train$unique_patient %in% id_thr),]
    bootstrap_slice <- as.numeric(as.matrix(data_thr$unique_patient))
    qmmd <- quantile_MMDl_bootstrap(Xthr, n_test, B_iterations, alpha, bootstrap_slice = bootstrap_slice, parallel = F)
    print(paste("id is ", id, " end bootstrap for alpha..."))
    ######################################
    threshold <- qmmd$threshold
    sigma <- qmmd$sigma
    power <- rep(0, length(alpha))
    for (id_val in unique(data_val$unique_patient)){
      Xval <- as.matrix(data_val[data_val$unique_patient == id_val,c("CCR7","CD4","CD45RA","CD3","HLADR","CD38","CD8")])
      Xval <- Xval[sample(nrow(Xval), n_test, replace = (nrow(Xval) < n_test)),]
      power <- power + test_MMDl_bootstrap(Xthr, Xval, threshold, sigma)
    }
    power <- power/length(unique(data_val$unique_patient))
    names(power) <- 1:length(alpha)
    if (power[1] > 0.05){
      best_alpha <- 0
    }else{
      best_alpha <- alpha[max(which(power <= 0.05))]
    }
    
    ### Do the tests! ###
    Xtrain <- as.matrix(data_train[,c("CCR7","CD4","CD45RA","CD3","HLADR","CD38","CD8")])
    bootstrap_slice <- as.numeric(as.matrix(data_train$unique_patient))
    ################################ Start parallelization to compute threshold for T bootstrap############################
    qmmd <- quantile_MMDl_bootstrap(Xtrain, n_test, B_iterations, alpha = best_alpha, bootstrap_slice = bootstrap_slice, parallel = F)
    print(paste("id is ", id, " end bootstrap for threshold..."))
    ######################################
    threshold <- qmmd$threshold
    sigma <- qmmd$sigma
    
    Xtest_0 <- as.matrix(data_test[data_test$Y == 0,c("CCR7","CD4","CD45RA","CD3","HLADR","CD38","CD8")])
    Xtest_1 <- as.matrix(data_test[data_test$Y == 1,c("CCR7","CD4","CD45RA","CD3","HLADR","CD38","CD8")])
    
    res_id <- data.frame(matrix(nrow = 0, ncol = 5))
    names(res_id) <- c("id", "theta",  "cell", "alpha", "power")
    for (theta in theta_choice){
      reject <- rep(0, length(alpha))
      for (r in 1:rep){
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
      res_id <- rbind(res_id, 
                      data.frame(id = id, theta = theta, cell = chosen_cell, alpha = best_alpha, power = reject/rep))
    }
    duration <- Sys.time() - start_time
    print(paste("cell type is :", chosen_cell, ", id is :", id, " best alpha is ", best_alpha, "computation has taken ", duration))
    print("")
    print("")
    saveRDS(res_index, paste("results/HIPC_MMD_maxsize/", chosen_cell, "_", id, ".rds", sep = ""))
    return(res_id)
  }
  saveRDS(res_cell, paste("results/ExpHIPC/mmd_", chosen_cell, ".rds", sep = ""))
}

parallel::stopCluster(cl = my.cluster)