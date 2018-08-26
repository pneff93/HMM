#' Normalized Backward Probability function
#'
#' @description This function calculates the normalized backward probability of the HMM states. It is important
#' to recall, that the provided weights come from the calculation of the forward probabilities one step ahead and
#' have to be entered as factor c 
#' In the general literatur also revered to as beta matrix.
#'
#' @param m number of likelihoods
#' @param N length of the supplied dataset
#' @param Gamma Gamma matrix
#' @param p vector of likelihood probabilities of the dataset
#' @param set index vector to align the p vector
#' @param c weights provided by the alpha calculation
#'
#' @return Returns the beta as a matrix with the dimensions m as rows
#' and the datapoints t as columns.



beta_function<-function( m, N, Gamma, p, set,c ){

  beta<-matrix(, ncol = N, nrow = m)
  beta[,N]<-rep(1, times=m)/c[N]

  for (t in (N-1):1){beta[,t] <- Gamma%*%diag(c(p[set+t]))%*%beta[,t+1]/c[t]}

  return(beta)
}


