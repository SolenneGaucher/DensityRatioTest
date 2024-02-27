### User defined splitting criteria for rpart
library(rpart)
source("utils/aux.R")

### Implements the functions init, eval, split to pass as arguments in rpart
### The objective is to find the partition that greedily maximizes sigma

### Initialisation function ### 
init_dr <- function(y, offset, parms, wt) {
  sfun <- function(yval, dev, wt, ylevel, digits) {
    paste("mean=", format(signif(yval, digits)),
          ",MSE=" , format(signif(dev/wt, digits)),
          sep = '')
  }
  environment(sfun) <- .GlobalEnv
  list(y = c(y), parms = parms, numresp = 1, numy = 1, summary = sfun)
}

### Evaluation function ###
eval_dr <- function(y, wt, parms){
  sigma <- compute_sigma((length(y)-sum(y)), sum(y), parms$n0, parms$n1, parms$epsilon0, parms$epsilon1)
  wmean <- sum(y)/length(y) 
  rss <- sum(wt*(y-wmean)^2) # necessary for computing the tree using rpart, but any criteria strictly decreasing will do
  list(label = sigma, deviance = rss)
}

################################### Splitting function ############################ 
split_dr <- function(y, wt, x, parms, continuous){
  if (! continuous){
    stop("Catagorical variables are not allowed")
  }
  n <- length(y)
  # number of points on the left and right sides of x[i] under f0 and f1
  N0_left <- cumsum((rep(1, n) - y)[-n])
  N0_right <- rev(cumsum((rep(1, n) - rev(y))[-n]))
  N1_left <- cumsum(y)[-n]
  N1_right <- rev(cumsum(rev(y))[-n])
  goodness <- sapply(1:(n-1), function(i) compute_delta_sigma(N0_left = N0_left[i], N0_right =  N0_right[i],
                                                              N1_left = N1_left[i], N1_right = N1_right[i], 
                                                              n0 = parms$n0, n1 = parms$n1, epsilon0 = parms$epsilon0, epsilon1 = parms$epsilon1))
  direction <- sapply(1:(n-1), function(i) sign(
    compute_sigma(N0_right[i], N1_right[i], n0 = parms$n0, n1 = parms$n1, epsilon0 = parms$epsilon0, epsilon1 = parms$epsilon1) - 
      compute_sigma(N0_left[i], N1_left[i], n0 = parms$n0, n1 = parms$n1, epsilon0 = parms$epsilon0, epsilon1 = parms$epsilon1)))
  return(list(goodness = goodness, direction =  direction)) # we set the most informative sub-bin as the left child
}

###### Group these functions as a method ################
method_dr <- function(){
  return(list(eval = eval_dr, split = split_dr, init = init_dr))
}