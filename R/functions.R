#' Fitting a Hidden Markov Model
#' @description Estimation of the transition probabilites, the initial state probabilites and the hidden state parameters of a Hidden Markov Model
#' by using the Direct Maximisation of the Likelihood or the EM-Algorithm
#' @param x, a sample of a Mixed Model
#' @param m, the number of states
#' @param method, choose between two different methods: "DM" as default, alternative "EM"
#' @param L1, likelihood of the first hidden state
#' @param L2, likelihood of the second hidden state
#' @param L3-L5, optional. likelihoods of the third, 4th and 5th hidden state
#' @param iterations, optional. Number of iterations for the EM-Algorithm
#' @param delta, optional. Stop criterion for the EM-Algorithm
#' @return The estimated parameters are rounded by 3 decimals and returned in a list.
#' @details 
#' This package is designed to estimate the hidden states of a HMM-Model, given the underlyning Likelihoods
#' of each state. It is important to support at least two Likelyhoods (L1,L2) for the function, which both
#' depend on a unknown parameter theta. 
#' 
#' In the method parameter the underlying estimation function is selected. With DM the HMM-function will
#' estimate the parameters via a direct maximisation of the given Likelihoods.
#' If EM is selected the HMM-function will use a Baum-Welch estimation algorithm to to compute the different states and the
#' estimation of the underlying parameters. 
#'
#' For more detailed explenation we recomend the source Hidden Markov Models for Times Series 
#' by Walter Zucchini, Iain MacDonald & Roland Langrock
#' 
#' 
#' The underlying functions are the HMM2 for the EM-algorithm and the HMM3 for the direct Maximisation
#' 
#' @export
#'
HMM<-function(x, m,method="DM", L1, L2, L3=NULL, L4=NULL, L5=NULL, iterations=NULL, delta=NULL){

#This is the head function, which seperates between the methods: 
  if (method == "DM"){
    output <- HMM3(x=x, m=m, L1=L1, L2=L2,L3=L3,L4=L4,L5=L5)
    return(output)
  } else if (method == "EM"){
    output <- HMM2(x=x, m=m, L1=L1, L2=L2, L3=L3, L4=L4, L5=L5, iterations=iterations, delta=delta)
    return(output)
    }else {
    warning("The supplied method is not available in the HMM function. Chooose between \"DM\" or \"EM\".", call. = T)
  }
}
