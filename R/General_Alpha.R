#' Normalized Forward Probability function
#'
#' @description This function calculates the normalized forward probability of 
#' the HMM states. The probabilities are implemented by calculating first the 
#' forward probabilities of each step and then they are weighted such that the 
#' sum of probabilities are equal to one. This procedure is done to prevent the 
#' threat of underflow. In the general literatur the forward probabilities are
#' also revered to as alpha matrix.
#'
#' @param m number of likelihoods
#' @param N length of the supplied dataset
#' @param delta Delta vector
#' @param Gamma Gamma matrix
#' @param p vector of likelihood probabilities of the dataset
#' @param set index vector to align the p vector
#'
#' @return Returns a list of parameters. The first list are the normalised 
#' alphas as a matrix with the dimension m as columns and the datapoints t as 
#' rows. The second part of the list are the corresponding weights.




alpha_function<-function( m, N, delta, Gamma, p, set ){

  alpha<-matrix(, ncol = m, nrow = N)
  weight <- c()
  weight[1] <- sum(t(delta)%*%diag(c(p[set])))
  alpha[1, ]<-(t(delta)%*%diag(c(p[set])))/weight[1]

  for (t in 2:N){
    prod <- alpha[t-1, ]%*%Gamma%*%diag(c(p[set+t-1]))
    weight[t] <- sum(prod)
    alpha [t, ] <- prod/weight[t]
  }
  out <- list(alpha, weight)
  return(out)
}
