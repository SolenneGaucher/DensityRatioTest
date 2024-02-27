rm(list=ls())
set.seed(0)
setwd("/workdir2/solenne.gaucher/FlowCytometry")
#setwd("/workdir2/solenne.gaucher/FlowCytometry")

# Load packages
library(doParallel)
library(partykit)
library(foreach)
library(lubridate)    
source("utils/user_defined_split.R")
source("utils/aux.R")
source("utils/test_dr.R")
source("utils/part_choice.R")


################################ Parameters ############################
n_sample_part <- 12 # 12 samples are used for choosing the partition
cell_choice <- c(2, 6, 8, 9)
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

############################################################
print(paste("We begin the experiment, time is ", now()))
for (chosen_cell in cell_choice){
  # label Y = 1 if cell == chosen_cell, 0 otherwise
  data_chosen_cell <- data_df
  data_chosen_cell$cell_type <- 1*(data_chosen_cell$cell_type == chosen_cell)
  names(data_chosen_cell)[names(data_chosen_cell) == 'cell_type'] <- "Y"
  
  n_test <- 10000
  id_part <- sample(unique(data_chosen_cell$unique_patient), n_sample_part)
  data_part <- data_chosen_cell
  n1_part <- sum(data_part$Y)
  n0_part <- nrow(data_part) - n1_part
    
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
  
  h <- compute_h_partition(partition, data_part, parms_part)
  sigma <- sum((h$h1k/h$h0k - 1)**2 * h0k)
  res <- data.frame(sigma = sigma)
  print(paste('cell is ', chosen_cell, ', and sigma is ', sigma))
  saveRDS(res, paste("results/HIPC_BDRT/estimated_signal", chosen_cell, "tot.rds", sep = ""))
}

