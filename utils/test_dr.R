### Auxilliary functions for the simulations
source("utils/aux.R")
library(partykit)

######################### Test T ###########################
test_dr <- function(data_train, data_test, parms, partition){
  # Returns the decision of the Density Ratio Test
  # data_test is the test dataset (type = data.frame)
  # data_train is the training dataset used to estimate h0K and h1k (type = data.frame)
  # partition is the partition (type = party)
  # parms are the corresponding parameters, and include:
  #   - n0: number of points from f0 in data_train
  #   - n1: number of points from f1 in data_train
  #   - epsilon0: threshold
  #   - epsilon1: threshold
  #   - u
  #   - t
  #   - K: size of the partition
  
  stats_train <- compute_h_partition(partition, data_train, parms)
  h0k <- stats_train$h0k
  h1k <- stats_train$h1k
  Omega <- stats_train$Omega
  rk <- h1k/h0k
  sigma <- sum((rk-1)**2 * h0k)
  n_test <- nrow(data_test)
  leaf_names <- unique(partition$fitted$`(fitted)`)
  Nk <- table(as.factor(predict(partition, newdata = data_test, type = "node")))
  for (node in leaf_names){
    if (!node %in% names(Nk)){
      Nk <- c(Nk, 0)
      names(Nk)[length(Nk)] <- node
    }
  }
  Nk <- Nk[order(names(Nk))]
  reject <- 0
  if (sum(Nk * (rk - 1))/n_test - 
      (sqrt(sigma) * (sqrt(10*parms$u*parms$K/parms$n0) + sqrt(6*parms$t/n_test)) + parms$t/(3*n_test) * max(abs(rk - 1)) + 3*parms$epsilon1*Omega) > 0){
    reject <- 1
  }
  return(reject)
}

######################### Approximate the distribution of S(X) under H0 by bootstrap ###########################
quantile_dr_bootstrap <- function(data, partition, parms, B_iterations = 1000, alpha = 0.05, bootstrap_slice = NULL){
  # Estimates the quantile of level alpha of the  statistic S(X) under H0
  # data is the dataset (type = data.frame)
  # partition is the partition (type = party)
  # parms are the corresponding parameters, and include:
  #   - n0: number of points from f0 for training
  #   - n1: number of points from f1 for training
  #   - n_test: number of points in the test sample
  #   - epsilon0: threshold
  #   - epsilon1: threshold
  #   - K: size of the partition
  
  T_distribution <- c()
  leaf_names <- unique(partition$fitted$`(fitted)`)
  n0 <- parms$n0
  n1 <- parms$n1
  n_test <- parms$n_test
  for (b in 1:B_iterations){
    # Create train and test dataset
    if (missing(bootstrap_slice)){
      bootstrap_test_index <- sample(which(data$Y == 0), parms$n_test, replace = TRUE)
      bootstrap_train_index0 <- sample(setdiff(which(data$Y == 0), bootstrap_test_index), n0, replace = TRUE)
      bootstrap_train_index1 <- sample(which(data$Y == 1), n1, replace = TRUE)
    }else{
      id_boostrap <- sample(unique(bootstrap_slice), 1)
      bootstrap_test_index <- sample(which(data$Y == 0 & bootstrap_slice == id_boostrap), parms$n_test, replace = TRUE)
      bootstrap_train_index0 <- sample(which(data$Y == 0 & bootstrap_slice != id_boostrap), n0, replace = TRUE)
      bootstrap_train_index1 <- sample(which(data$Y == 1 & bootstrap_slice != id_boostrap), n1, replace = TRUE)
    }
    bootstrap_train0 <- data[bootstrap_train_index0,]
    bootstrap_train1 <- data[bootstrap_train_index1,]
    bootstrap_test <- data[bootstrap_test_index,]
    
    # Compute the estimate d parameters and signal for the test with bootstrap samples
    stats_bootstrap <- compute_h_partition(partition, rbind(bootstrap_train0, bootstrap_train1), parms)
    h0k <- stats_bootstrap$h0k
    h1k <- stats_bootstrap$h1k
    rk <- h1k/h0k
    Nk <- table(as.factor(predict(partition, newdata = bootstrap_test, type = "node")))/parms$n_test
    for (node in leaf_names){
      if (!node %in% names(Nk)){
        Nk <- c(Nk, 0)
        names(Nk)[length(Nk)] <- node
      }
    }
    Nk <- Nk[order(names(Nk))]
    T_distribution <- c(T_distribution, mean((rk-1)*Nk))
  }
  return(quantile(T_distribution, 1-alpha))
}

test_dr_bootstrap <- function(data_train, data_test, partition, parms, threshold = NULL, alpha = 0.05, B_iterations = 1000, bootstrap_slice = NULL){
  # Runs the Boostraped Density Ratio Test
  # data is the dataset (type = data.frame)
  # partition is the partition (type = party)
  # parms are the corresponding parameters, and include:
  #   - n0: number of points from f0 for training
  #   - n1: number of points from f1 for training
  #   - epsilon0: threshold
  #   - epsilon1: threshold
  #   - K: size of the partition
  
  if (missing(threshold)){
    threshold <- quantile_dr_bootstrap(data_train, partition, parms, B_iterations, alpha, bootstrap_slice)
  }
  
  stats_train <- compute_h_partition(partition, data_train, parms)
  leaf_names <- unique(partition$fitted$`(fitted)`)
  h0k <- stats_train$h0k
  h1k <- stats_train$h1k
  rk <- h1k/h0k
  n_test <- dim(data_test)[1]
  Nk <- table(as.factor(predict(partition, newdata = data_test, type = "node")))/n_test
  for (node in leaf_names){
    if (!node %in% names(Nk)){
      Nk <- c(Nk, 0)
      names(Nk)[length(Nk)] <- node
    }
  }
  Nk <- Nk[order(names(Nk))]
  return(1*(mean((rk-1)*Nk) > threshold))
}
