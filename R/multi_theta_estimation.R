#' Multi Theta Estimation function
#'
#' @description This function calculates part of the global log-likelihood that 
#' is only dependent on the Theta value. Due to its proportionality, it is 
#' therefore optimal for the maximisation of the Theta values and will be used 
#' by the EM-algorithm. For the multi_HMM_EM(), small changes where made to 
#' calculate the index of the Theta alues right.
#' 
#' @param m number of likelihoods
#' @param N length of the supplied dataset
#' @param u output matrix of the u_function
#' @param x a sample of a Hidden Markov Model 
#' @param theta Theta vector
#' @param L1 likelihood of the first hidden state
#' @param L2 likelihood of the second hidden state
#' @param L3 optional. likelihoods of the third hidden state
#' @param L4 optional. likelihoods of the 4th hidden state
#' @param L5 optional. likelihoods of the 5th hidden state
#' @param start_index index paramter to assign the right amount of thetas 
#' to the likelihoods 
#'
#'@details For more detailed explanation we recommend the source Hidden Markov 
#'Models for Times Series by Walter Zucchini, Iain MacDonald & Roland Langrock, 
#'especially page 72.
#'
#'
#' @return returns the corresponding likelihood
#'


multi_theta_estimation <- function( m, N, u, x, theta, L1, L2, L3 = NULL, 
                                  L4 = NULL,L5 = NULL, start_index ){
  
  #For each state, the individual likelihood is logarithmized and summed with 
  #the coresponding u-vector over all timepoints. The contributions of the 
  #states are then added together to a global likelihood that needs to be 
  #maximized 
  
  #Due to the fact that the maximization is dependent on the number of likelihoods
  #we query for the given number of states 
  
  #The index start_index is used to rightly identify the first and the last 
  #theta-value of the grouping vector Theta (see description in multi_HMM_EM())

  l3<-0
  

  
  #For the first likelihood
  l3 <- l3 + u[, 1] %*% (log(L1(x, theta[start_index[1]:(start_index[2]-1)])))
  #For the second likelihood
  l3 <- l3 +u[, 2] %*% (log(L2(x, theta[start_index[2]:(start_index[3]-1)])))
  
  #If we have more than two likelihoods we can use this algorithm
  if (m>2){
    l3 <- l3 +u[, 3 ]%*% (log(L3(x, theta[start_index[3]:(start_index[4]-1)])))
  }
  if (m>3) {
    l3 <- l3 +u[, 4] %*% (log(L4(x, theta[start_index[4]:(start_index[5]-1)])))
  }
  if (m>4){
    l3 <- l3 +u[, 5] %*% (log(L5(x, theta[start_index[5]:(start_index[6]-1)])))
  }
  return(l3*-1)
}

