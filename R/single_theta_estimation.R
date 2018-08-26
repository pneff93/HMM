#' Theta Estimation function
#'
#' @description This function calculates part of the global log-likelihood that is only dependend 
#' on the theta value. Due to its proportionality, it is therefore optimal for the maximisation 
#' of the theta-values and will be used by the EM-algorithm. 
#' 
#' 
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
#' @details For more detailed explanation we recommend the source Hidden Markov Models for Times Series
#' by Walter Zucchini, Iain MacDonald & Roland Langrock, especially page 72.
#'
#'
#' @return returns the corresponding likelihood
#'


theta_estimation<-function( m, N, u, x, theta, L1, L2, L3=NULL, L4=NULL, L5=NULL ){


  
  #For each state, the individual likelihood is logarithmized an summed with 
  #the coresponding u-vector over all timepoints. 
  #The contributions of the states are then added together to a global likelihood that needs to be maximized 
  
  #Due to the fact that the maximization is dependen on the number of likelihoods
  #we query for the given number of states 
  
  l3<-0
  
  #For the first likelihood
  l3 <- l3 + u[,1]%*%(log(L1(x,theta[1])))
  #For the second likelihood
  l3 <- l3 +u[,2]%*%(log(L2(x,theta[2])))

  #If we have more than two likelihoods we add their their contributions 
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

