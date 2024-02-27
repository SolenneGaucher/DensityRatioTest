### Auxilliary functions for the simulations
library(kernlab)

MMDl <- function(x, y, sigma = 0.5){
  # Computes the MMDl statistic
  # x and y are two matrices 
  # Sample for X
  l <- floor(min(nrow(x), nrow(y))/2)
  idx <- sample(nrow(x), 2 * l)
  x1 <- x[idx[seq_len(l)], ]
  x1 <- as.matrix(x1)
  norms_x1 <- rowSums(x1^2)
  x2 <- x[idx[seq(l + 1, 2 * l)], ]
  x2 <- as.matrix(x2)
  norms_x2 <- rowSums(x2^2)
  
  # Sample for Y
  idy <- sample(nrow(y), 2 * l)
  y1 <- y[idy[seq_len(l)], ]
  y1 <- as.matrix(y1)
  norms_y1 <- rowSums(y1^2)
  y2 <- y[idx[seq(l + 1, 2 * l)], ]
  y2 <- as.matrix(y2)
  norms_y2 <- rowSums(y2^2)
  
  # Kernels
  rbf <- rbfdot(sigma = sigma)
  Kx <- as.vector(kernlab::kernelFast(rbf, x1, x2, norms_x1))
  Ky <- as.vector(kernlab::kernelFast(rbf, y1, y2, norms_y1))
  Kxy <- as.vector(kernlab::kernelFast(rbf, x1, y1, norms_x1))
  Kyx <- as.vector(kernlab::kernelFast(rbf, y2, x2, norms_y2)) 
  
  # The MMD^2_u linear statistic.
  term1 <- (1 / l) * sum(Kx)
  term2 <- (1 / l) * sum(Ky)
  term3 <- (1 / l) * sum(Kxy)
  term4 <- (1 / l) * sum(Kyx)
  return(term1 + term2 - term3 - term4)
}

quantile_MMDl_bootstrap <- function(Xtrain, n_test, B_iterations = 1000, 
                                    alpha = 0.05, bootstrap_slice = NULL, parallel = FALSE){
  # Computes the MMDl quantiles of level alpha
  # Xtrain : data matrix
  # n_test : size of the test sample, i.e. of the samples to use in the bootstrap
  T_distribution <- c()
  n_train <- nrow(Xtrain)
  # Sample for estimating the bandwidth = median distance
  sigma <- median(dist(Xtrain[sample(1:nrow(Xtrain), 1000, replace = (nrow(Xtrain) < 1000)), ]))
  if (!parallel){
    for (b in 1:B_iterations){
      if (missing(bootstrap_slice)){
        bootstrap_test_index <- sample(1:n_train, n_test, replace = FALSE)
        bootstrap_train_index <- sample((1:n_train)[-bootstrap_test_index], n_test, replace = (n_train < 2*n_test))
      }else{
        id_boostrap <- sample(unique(bootstrap_slice), 1)
        bootstrap_test_index <- sample(which(bootstrap_slice == id_boostrap), n_test, replace = T)
        bootstrap_train_index <- sample(which(bootstrap_slice != id_boostrap), n_test, replace = T)
      }
      bootstrap_test <- Xtrain[bootstrap_test_index,]
      bootstrap_train <- Xtrain[bootstrap_train_index,]
      T_distribution <- c(T_distribution, MMDl(bootstrap_train, bootstrap_test, sigma))
    }
  }else{
    T_distribution <- foreach (b = 1:B_iterations, .combine = rbind) %dopar% {
      source("test_MMD.R")
      if (missing(bootstrap_slice)){
        bootstrap_test_index <- sample(1:n_train, n_test, replace = FALSE)
        bootstrap_train_index <- sample((1:n_train)[-bootstrap_test_index], n_test, replace = (n_train < 2*n_test))
      }else{
        id_boostrap <- sample(unique(bootstrap_slice), 1)
        bootstrap_test_index <- sample(which(bootstrap_slice == id_boostrap), n_test, replace = T)
        bootstrap_train_index <- sample(which(bootstrap_slice != id_boostrap), n_test, replace = T)
      }
      bootstrap_test <- Xtrain[bootstrap_test_index,]
      bootstrap_train <- Xtrain[bootstrap_train_index,]
      return(MMDl(bootstrap_train, bootstrap_test, sigma))
    }
  }
  return(list(threshold = quantile(T_distribution, 1-alpha), sigma = sigma))
}

test_MMDl_bootstrap <- function(data_train, data_test, threshold, sigma){
  n_train <- nrow(data_train)
  n_test <- nrow(data_test)
  data_train_MMDl <- data_train[sample(1:n_train, n_test, replace = FALSE),]
  return(1*(MMDl(data_train_MMDl, data_test, sigma) > threshold))
}