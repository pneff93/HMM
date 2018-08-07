#' Normalized Backward Probability function
#' 
#' @description This function calculates the normalized backward probability of the HMM states. It is important
#' to recall, that the provided weights come from the calculation of the forward probabilities one step ahead
#' In the general literatur also revered to as beta matrix
#' 
#' @param m number of Likelihoods
#' @param N length of the supplied dataset
#' @param Gamma Gamma matrix 
#' @param p vector of Likelihood probabilities of dataset
#' @param set index vector to align the p vector 
#' @param c weights provided by the alpha calculation
#' 
#' @return Returns the beta as a matrix with the dimensions m rows
#' and the datapoints t as collumns



beta_function<-function( m, N, Gamma, p, set,c ){
  
  beta<-matrix(, ncol = N, nrow = m)
  beta[,N]<-rep(1, times=m)/c[N]
  
  for (t in (N-1):1){beta[,t] <- Gamma%*%diag(c(p[set+t]))%*%beta[,t+1]/c[t]}
  
  return(beta)
}


