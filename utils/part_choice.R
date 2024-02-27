### Auxilliary functions for the simulations
source("utils/aux.R")
library(rpart)
library(partykit)
library(kernlab)
library(pbapply)

######################### Construct a sequence of partitions, choose the best one ###########################
# PruneTree
prune_dr <- function(full_tree, data, parms, Kmin = 3, verbose = F){
  # Returns the best partition among the partition sequence described by full_tree
  # data is the dataset used for choosing the best partition in the sequence
  # full¨tree is the partition sequence (type = rpart)
  # parms are the corresponding parameters, and include:
  #   - n0: number of points from f0 in data
  #   - n1: number of points from f1 in data
  #   - epsilon0: threshold
  #   - epsilon1: threshold
  
  # create copies of objects
  snip_tree <- full_tree
  snip_data <- data
  
  # store intermediate results
  snip_order <- c()
  snip_sigma <- c()
  partition_size <- c()
  count <- 1
  
  if (nrow(snip_tree$frame) <= Kmin){
    data$leaf <- as.factor(row.names(snip_tree$frame)[snip_tree$where])
    snip_data$leaf <- as.factor(row.names(snip_tree$frame)[snip_tree$where])
    Ntable <- table(snip_data$Y, snip_data$leaf) # count the number of points in each leave for X0part and X1part
    snip_sigma <- c(sum(sapply(1:ncol(Ntable), function(k) compute_sigma(Ntable[1,k], Ntable[2,k], 
                                                                         n0 =parms$n0, n1 =parms$n1, 
                                                                         epsilon0 =parms$epsilon0, epsilon1 =parms$epsilon1))), snip_sigma)
    return(list(sigma = sigma, pruned_tree = full_tree, K = nrow(full_tree$frame), partition_size = c(nrow(full_tree$frame))))
  }
  
  # compute the order for snipping the nodes
  while (nrow(snip_tree$frame) > Kmin){ # there are still more than Kmin leafs/bins in the partition
    
    # Compute the signal sigma corresponding to the partition snip_tree
    snip_data$leaf <- as.factor(row.names(snip_tree$frame)[snip_tree$where])
    Ntable <- table(snip_data$Y, snip_data$leaf) # count the number of points in each leave for X0part and X1part
    snip_sigma <- c(sum(sapply(1:ncol(Ntable), function(k) compute_sigma(Ntable[1,k], Ntable[2,k], 
                                                                    n0 =parms$n0, n1 =parms$n1, 
                                                                    epsilon0 =parms$epsilon0, epsilon1 =parms$epsilon1))), snip_sigma)
    partition_size <- c(ncol(Ntable), partition_size)
    # Find the least informative pair of adjacent terminal nodes
    leaf_pairs <- c()
    delta_sigma <- c()
    for (i in 1:(nrow(snip_tree$frame)-1)){
      if (snip_tree$frame$var[i] == "<leaf>" && snip_tree$frame$var[i+1] == "<leaf>" &&
          as.numeric(row.names(snip_tree$frame)[i]) == as.numeric(row.names(snip_tree$frame)[i+1])-1){ # nodes correspond to adjacent leafs
        leaf_pairs <- c(leaf_pairs, as.numeric(row.names(snip_tree$frame)[i]))
        delta_sigma <- c(delta_sigma, compute_delta_sigma(Ntable[1,row.names(snip_tree$frame)[i]], Ntable[1,as.character(as.numeric(row.names(snip_tree$frame)[i])+1)],
                                                          Ntable[2,row.names(snip_tree$frame)[i]], Ntable[2,as.character(as.numeric(row.names(snip_tree$frame)[i])+1)],
                                                         parms$n0,parms$n1,parms$epsilon0,parms$epsilon1))
      }
    }
    snip_node <- leaf_pairs[which.min(delta_sigma)]/2 # we snip the tree under the parent node of the less informative node pair
    snip_order <- c(snip_order, snip_node)
    snip_tree <- snip.rpart(snip_tree, snip_node)
    if (verbose){
      print(paste(count, " th aggregation : node number ", snip_node, "is sniped, decrease in signal : ", min(delta_sigma)))
    }
    count <- count + 1
  }
  
  # Choose the tree size
  K <- partition_size[max(which(snip_sigma/partition_size >= 0.9*max(snip_sigma/partition_size)))]
  sigma <- snip_sigma[max(which(snip_sigma/partition_size >= 0.9*max(snip_sigma/partition_size)))]
  
  # Prune the tree to obtain for the desired number of leaves K
  pruned_tree <- full_tree
  count <- 1
  while (sum(pruned_tree$frame$var == "<leaf>") > K){
    pruned_tree <- snip.rpart(pruned_tree, snip_order[count])
    count <- count + 1
  }
  return(list(sigma = sigma, pruned_tree = pruned_tree, K = K, snip_sigma = snip_sigma, partition_size = partition_size))
}
