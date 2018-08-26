#' Likelihood of the Hidden Markov Model
#'
#' @description This function calculates the log-likelihood of the HMM mode.
#' For this, the scaled forward probabilities are computed. In respect of the LH
#' function, the multiLH function incorporates the multi Thetas for the 
#' indivudal likelihoods 
#' 
#' @param factor Input of variables that are unrestricted
#' @param x a sample of a Hidden Markov Model
#' @param m the number of states
#' @param L1 likelihood of the first hidden state
#' @param L2 likelihood of the second hidden state
#' @param L3 optional. likelihood of the third hidden state
#' @param L4 optional. likelihood of the 4th hidden state
#' @param L5 optional. likelihood of the 5th hidden state
#' @param start_index index parameter to assign the Thetas to their corresponding 
#' individual hidden state likelihood.
#'
#' @return negative Likelihood
#' @details This function computes the log-likelihood of the forward 
#' probabilities of the HMM. Given the fact that the inputed factor vector is 
#' not restricted, we need to apply a transformation to transform the factor
#' variables to our suitable canditates Delta,Gamma and Theta. For this we apply 
#' the function trans().
#' 
#' The likelihood is constructed using the scaling/normalizing of the forward 
#' probabilities such that for each time point the sum of the foward probabily is
#' equal to one. The scaling is neccesary to tackle the underflow problem, that 
#' arrises with an increasing sample size.
#' 
#' The multiLH function only differes in the computation of the likelihood 
#' probability vector p, because we imput multiple Thetas in the corresponding 
#' individual likelihoods. 
#' 
#'



multiLH <- function (factor, x, m, L1, L2, L3 = NULL, L4 = NULL, L5 = NULL,
                     start_index){
  
  ##########################
  #Transformation
  
  #We first have to transform the factors without constrains into our 
  #Delta/Gamma/Theta values with constrains
  out <- trans(factor, m)
  delta <- out[[1]]
  gamma <- out[[2]]
  theta <- out[[3]]
  
  N <- length(x)
  
  
  #likelihoods
  
  #We combine the individual likelihoods for each timepoint to on general 
  #probability vector. The index set defines every cut, where a new likelihood
  #beginns. This measur increase the performance of our calculations drastically
  #(instead of calculating a huge diagonal matrix)
  
  p1<-L1(x, theta[start_index[1]:(start_index[2]-1)])
  p2<-L2(x, theta[start_index[2]:(start_index[3]-1)])
  p<-c(p1, p2)
  
  if(!is.null(L3)){
    p3<-L3(x, theta[start_index[3]:(start_index[4]-1)])
    p<-c(p, p3)
  }
  if(!is.null(L4)){
    p4<-L4(x, theta[start_index[4]:(start_index[5]-1)])
    p<-c(p, p4)
  }
  if(!is.null(L5)){
    p5<-L5(x, theta[start_index[5]:(start_index[6]-1)])
    p<-c(p, p5)
  }    
  
  set<-seq(1, length(p)-N+1, length.out = m)
  
  ##########################
  #Computation of normalized forward probabilities
  
  #Computation of the log-likelihood with normalized alphas to tackle the 
  #underflow problem
  
  nalpha<- matrix(, nrow = m, ncol = N)
  v <- delta%*%diag(c(p[set]))
  u <- sum(v)
  l <- log(u)
  nalpha [, 1] <- t(v/u)
  
  
  for (t in 2:N){
    v <- nalpha[, t-1]%*%gamma%*%diag(c(p[set+t-1]))
    u <- sum(v)
    l <- l + log(u)
    nalpha[, t]<- t(v/u)
    
  }
  return (-1*l)
}