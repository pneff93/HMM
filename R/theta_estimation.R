#' Theta Estimation function
#'
#' @description This function calculates the third part of the overall likelihood function. 
#' This likelihood part focusses on the hidden states parameters and will be optimized later on.
#' It is estimated as a sum over all time periods and states of the vector u times the corresponding log-likelihood.
#'
#' @param m number of likelihoods
#' @param N length of the supplied dataset
#' @param u output matrix of the u_function
#' @param x a sample of a Mixed Model
#' @param theta theta vector
#' @param L1 likelihood of the first hidden state
#' @param L2 likelihood of the second hidden state
#' @param L3 optional. likelihoods of the third hidden state
#' @param L4 optional. likelihoods of the 4th hidden state
#' @param L5 optional. likelihoods of the 5th hidden state
#'
#' @return returns the corresponding likelihood
#'


theta_estimation<-function( m, N, u, x, theta, L1, L2, L3=NULL, L4=NULL, L5=NULL ){

  l3<-0

  #For the first likelihood
  l3 <- l3 + u[,1]%*%(log(L1(x,theta[1])))
  #For the second likelihood
  l3 <- l3 +u[,2]%*%(log(L2(x,theta[2])))

  #If we have more than two likelihoods we can use this algorithm
  if (m>2){
    l3 <- l3 +u[,3]%*%(log(L3(x,theta[3])))
  }
  if (m>3) {
    l3 <- l3 +u[,4]%*%(log(L4(x,theta[4])))
  }
  if (m>4){
    l3 <- l3 +u[,5]%*%(log(L5(x,theta[5])))
  }
  return(l3*-1)
}

