#' Likelihood of the Hidden Markov Model
#'
#' @description This function calculates the log-likelihood of the HMM model by normalizing the forward-probabilities.
#'
#' @param factor Input of variables that are unrestricted
#' @param x a sample of a Mixed Model
#' @param m the number of states
#' @param L1 likelihood of the first hidden state
#' @param L2 likelihood of the second hidden state
#' @param L3 optional. likelihood of the third hidden state
#' @param L4 optional. likelihood of the 4th hidden state
#' @param L5 optional. likelihood of the 5th hidden state
#'
#' @return negative Likelihood
#' @details This function computes the log-likelihood of the forward probabilities of the HMM given the factor values
#' Due to the fact that the Sigma and Gamma matrices are in a form restricted, the transformation function is used to 
#' transform the factor variables in suitable canidates. 
#' 
#' The likelihood is constructed with normalizing the alpha vectors.
#' 
#'






LH <- function (factor,x,m, L1, L2,L3=NULL,L4=NULL,L5=NULL){
  
  
  #We first have to transform the factors without constrains into our Sigma/Gamma/Theta
  #values with constrains
  out <- trans(factor,m)
  sigma <- out[[1]]
  gamma <- out[[2]]
  theta <- out[[3]]
  
  T <- length(x)
  
  
  #likelihoods
  p1<-L1(x, theta[1])
  p2<-L2(x, theta[2])
  p<-c(p1,p2)
  
  if(!is.null(L3)){
    p3<-L3(x, theta[3])
    p<-c(p,p3)
  }
  if(!is.null(L4)){
    p4<-L4(x, theta[4])
    p<-c(p,p4)
  }
  if(!is.null(L5)){
    p5<-L5(x, theta[5])
    p<-c(p,p5)
  }
  
  set<-seq(1, length(p)-T+1, length.out = m)
  
  #Computation of the log-likelihood with normalized alphas to tackle the underflow problem
  
  #normalized alpha
  nalpha<- matrix(,nrow=m, ncol=T)
  v <- sigma%*%diag(c(p[set]))
  u <- sum(v)
  l <- log(u)
  nalpha [,1] <- t(v/u)
  
  
  for (t in 2:T){
    v <- nalpha[,t-1]%*%gamma%*%diag(c(p[set+t-1]))
    u <- sum(v)
    l <- l + log(u)
    nalpha[,t]<- t(v/u)
    
  }
  return (-1*l)
}