### Auxilliary functions for the simulations
library(rpart)
library(partykit)
library(kernlab)

### Implements the functions init, eval, split to pass as arguments in rpart
### The objective is to find the partition that greedily maximizes sigma

### Simple auxilliary functions used to compute the thresholded density estimates ### 
compute_h0 <- function(h0, h1, epsilon0, epsilon1){
  # h0 and h1 are the (non-thresholded) frequencies
  # epsilon0 and epsilon1 are the thresholds
  
  rep <- h0
  if (h1 < epsilon1){
    if (h0 < epsilon1){
      rep <- min(3*epsilon1,1)
    }
  }else{
    if (h0 < epsilon0){
      rep <- min(3*epsilon0,1)
    }
  }
  return(rep)
}

compute_h1 <- function(h1, epsilon1){
  # h1 is the (non-thresholded) frequency
  # epsilon1 is the thresholds
  
  rep <- h1
  if (h1 < epsilon1){
    rep <-  min(3*epsilon1,1)
  }
  return(rep)
}

compute_sigma <- function(N0k, N1k, n0, n1, epsilon0, epsilon1){
  # Compute the signal in one bin
  # N0K and N1k are the number of points from f0 and f1 in the bin
  # n0 and n1 the total number of points from f0 and f1
  # epilon0 and epsilon1 are the thresholds
  
  h0 <- compute_h0(N0k/n0, N1k/n1, epsilon0, epsilon1)
  h1 <- compute_h1(N1k/n1, epsilon1)
  return((h1/h0-1)**2 * h0)
}

compute_delta_sigma <- function(N0_left, N0_right, N1_left, N1_right, n0, n1, epsilon0, epsilon1){
  # Compute the increase in signal for splitting a bin in two
  # N0K_left and N1k_left  are the number of points from f0 and f1 in the left bin
  # N0K_right and N1k_right  are the number of points from f0 and f1 in the rsight bin
  # n0 and n1 the total number of points from f0 and f1
  # epilon0 and epsilon1 are the thresholds
  
  sigma_left <- compute_sigma(N0_left, N1_left, n0, n1, epsilon0, epsilon1)
  sigma_right <- compute_sigma(N0_right, N1_right, n0, n1, epsilon0, epsilon1)
  sigma_tot <- compute_sigma(N0_left + N0_right, N1_left + N1_right, n0, n1, epsilon0, epsilon1)
  return(sigma_left + sigma_right - sigma_tot)
}


### Auxilliary functions used to estimates the statistics h and Omega ### 
compute_h_partition <- function(partition, data, parms){
  # Compute the density estimated probabilities h0k and h1k for the nodes 1, ..., K
  # partition is the partition (type = party)
  # data is the training dataset used to estimate h0K and h1k (type = data.frame)
  # parms are the corresponding parameters, and include:
  #   - n0: number of points from f0 in data
  #   - n1: number of points from f1 in data
  #   - epsilon0: threshold
  #   - epsilon1: threshold
  
  leaf_names <- unique(partition$fitted$`(fitted)`)
  data$leaf <- as.factor(predict(partition, newdata = data, type = "node"))
  Ntable_est <- table(data$Y, data$leaf) # count the number of points in each leave for X0part and X1part
  for (node in leaf_names){
    if (!node %in% colnames(Ntable_est)){
      newcol <- matrix(data = c(0,0), ncol = 1)
      colnames(newcol) <- node
      Ntable_est<- cbind(Ntable_est, newcol)
    }
  }
  Ntable_est <- Ntable_est[,order(colnames(Ntable_est))]
  h0k <- sapply(1:ncol(Ntable_est), function(k) compute_h0(Ntable_est[1,k]/parms$n0, Ntable_est[2,k]/parms$n1, 
                                                                epsilon0 = parms$epsilon0, epsilon1 = parms$epsilon1))
  names(h0k) <- colnames(Ntable_est)
  h1k <- sapply(1:ncol(Ntable_est), function(k) compute_h1(Ntable_est[2,k]/parms$n1, epsilon1 = parms$epsilon1))
  names(h1k) <- colnames(Ntable_est)
  Omega01 <- sum(Ntable_est[2,]/parms$n1 <= parms$epsilon1 & Ntable_est[1,]/parms$n0 <= parms$epsilon1)
  Omega1 <- sum(Ntable_est[2,]/parms$n1 <= parms$epsilon1 & Ntable_est[1,]/parms$n0 > parms$epsilon1)
  return(list(h0k = h0k, h1k = h1k, Omega = Omega01 + Omega1))
}
