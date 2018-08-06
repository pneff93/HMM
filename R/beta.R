#' Backward Probability function
#' 
#' @description This function calculates the backward probability of the HMM states. In the general 
#' literatur also revered to as beta matrix
#' 
#' @param m number of Likelihoods
#' @param N length of the supplied dataset
#' @param Gamma Gamma matrix 
#' @param p vector of Likelihood probabilities of dataset
#' @param set index vector to align the p vector 
#' 
#' @return Returns the beta as a matrix with the dimensions m rows
#' and the datapoints t as collumns



beta_function<-function( m, N, Gamma, p, set ){
  
  beta<-matrix(, ncol = N, nrow = m)
  beta[,N]<-rep(1, times=m)
  
  for (t in (N-1):1){
    beta[,t] <- Gamma%*%diag(c(p[set+t]))%*%beta[,t+1]
  }
  
  return(beta)
}
