#' Fitting a Hidden Markov Model
#' @description Estimation of the transition probabilites, the initial state probabilites and the hidden state parameters of a Hidden Markov Model
#' by using the Direct Maximisation of the likelihood or the EM-Algorithm.
#' @param x a sample of a Mixed Model
#' @param m the number of states
#' @param method choose between two different methods: "DM" as default, alternative "EM"
#' @param L1 likelihood of the first hidden state
#' @param L2 likelihood of the second hidden state
#' @param L3 optional. likelihood of the third hidden state
#' @param L4 optional. likelihood of the 4th hidden state
#' @param L5 optional. likelihood of the 5th hidden state
#' @param iterations optional. number of iterations for the EM-Algorithm
#' @param delta optional. stop criterion for the EM-Algorithm
#' @param decoding if parameter set TRUE the function returns the most propable paths via local and global decoding
#' 
#' 
#' @return The estimated parameters are rounded by 3 decimals and returned in a list
#' @details
#' This package is designed to estimate the hidden states of a HMM-Model, given the underlying likelihoods
#' of each state. It is important to support at least two likelihoods (L1, L2) for the function, which both
#' depend on one unknown parameter theta.
#'
#' In the method parameter the underlying estimation function is selected. With DM the HMM-function will
#' estimate the parameters via a direct maximisation of the given likelihoods.
#' If EM is selected the HMM-function will use a Baum-Welch estimation algorithm to compute the different states and the
#' estimation of the underlying parameters.
#'
#' For more detailed explanation we recommend the source Hidden Markov Models for Times Series
#' by Walter Zucchini, Iain MacDonald & Roland Langrock.
#'
#'
#' The underlying functions are the HMM2 for the EM-Algorithm and the HMM3 for the Direct Maximisation
#'
#' @export
#'
HMM<-function(x, m,method="DM", L1, L2, L3=NULL, L4=NULL, L5=NULL, iterations=NULL, delta=NULL,decoding=FALSE){

#This is the head function, which seperates between the methods:
  if (method == "DM"){
    output <- HMM3(x=x, m=m, L1=L1, L2=L2,L3=L3,L4=L4,L5=L5)
    
    #Decoding:
    if (decoding==TRUE){
      dec <-decode(x=x,m=m, L1=L1, L2=L2,L3=L3,L4=L4,L5=L5,gamma=output$Gamma, sigma=output$Sigma,theta=output$Theta,multi=FALSE)
      output <- (append(output,dec)) 
    }
    return(output)
  } else if (method == "EM"){
    output <- HMM2(x=x, m=m, L1=L1, L2=L2, L3=L3, L4=L4, L5=L5, iterations=iterations, delta=delta)
    
    #Decoding:
    if (decoding==TRUE){
      dec <-decode(x=x,m=m, L1=L1, L2=L2,L3=L3,L4=L4,L5=L5,gamma=output$Gamma, sigma=output$Sigma,theta=output$Theta,multi=FALSE)
      output <- append(output,dec)
    } 
    
    return(output)
    }else {
    warning("The supplied method is not available in the HMM function. Chooose between \"DM\" or \"EM\".", call. = T)
  }
}
