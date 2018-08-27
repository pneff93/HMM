#' Fitting a Hidden Markov Model to multifactor Likelihoods
#' 
#' @description Estimation of the transition probabilites, the initial state 
#' probabilites and the hidden state parameters of a Hidden Markov Model 
#' by using the Direct Maximisation of the likelihood or the Baum-Welch 
#' Algorithm
#' @param x sample of a Hidden Model 
#' @param theta list with initial number of Likelihood parameters (see details)
#' @param m the number of states
#' @param method choose between two different methods: "DM" as default, 
#' alternative "EM"
#' @param L1 likelihood of the first hidden state
#' @param L2 likelihood of the second hidden state
#' @param L3 optional. likelihood of the third hidden state
#' @param L4 optional. likelihood of the 4th hidden state
#' @param L5 optional. likelihood of the 5th hidden state
#' @param iterations optional. number of iterations for the EM-Algorithm
#' @param DELTA optional. stop criterion for the EM-Algorithm
#' @param decoding if parameter set TRUE the function returns the most 
#' probable paths via local and global decoding
#' 
#' @return Returns the Delta vector, Gamma matrix and the Thetas of the
#' Likelihoods rounded by three decimals. If "EM" is selected the function also
#' returns the number of iterations and the DELTA. 
#' 
#' @details
#' This package is designed to estimate the hidden states of a HMM-Model, given
#' the underlying likelihoods of each state. It is important to support at
#' least two likelihoods (L1, L2) for the function, which both depend on
#' multiple unkown thetas. See examples for a suitable structur of the
#' likelihood.
#' 
#' The multi_HMM() function is able to calculate with multiple Theta values
#' for the individual likelihoods. For each likelihood the right number initial
#' starting parameter has to be set, in order to compute the estimation of the 
#' corresponding Thetas.
#' For each Likelihood the starting values must be in the format of a vector, 
#' which is then saved as a list element. 
#' 
#' e.g.: theta[[i]] <- c(parameter1, parameter2, ...)
#' 
#' The function then extracts the right number of parameters per likelihood and 
#' optimizes the values. 
#' 
#' 
#'
#' Choose with "method" the underlying estimation function. If DM is selected, 
#' the HMM-function will estimate the parameters via a direct maximisation of 
#' the given likelihoods.
#' If EM is selected the HMM-function will use a Baum-Welch estimation algorithm
#' to compute the different states and the estimation of the underlying 
#' parameters.
#'
#' For more detailed explanation we recommend the source Hidden Markov Models 
#' for Times Series by Walter Zucchini, Iain MacDonald & Roland Langrock.
#'
#'
#' The underlying functions are the multi_HMM_EM for the EM-Algorithm and the 
#' multi_HMM_DM for the Direct Maximisation.
#'
#'@import stats
#'
#'@seealso For Hidden Markov Models with only one theta per likelihood, please
#'refer to \code{\link{single_HMM}} 
#'
#'
#' @export
#' 
#' 
#' 
#' @example R/examples/multi_Norm_threeStates.R
#' 
#' 
#'
multi_HMM <- function( x, theta, m, method = "DM", L1, L2, L3 = NULL,
                         L4 = NULL, L5 = NULL, iterations = NULL,
                         DELTA = NULL, decoding = FALSE ){
  
  #This is the head function, which seperates between the methods:
  if (method == "DM"){
    output <- multi_HMM_DM( x = x, theta = theta, m = m, L1 = L1, L2 = L2, L3 = L3,
                        L4 = L4, L5 = L5 )
    
    #Decoding
    if (decoding==TRUE){
      dec <- decode( x = x, m = m, L1 = L1, L2 = L2, L3 = L3, L4 = L4, L5 = L5,
                   gamma = output$Gamma, delta = output$Delta,
                   theta = output$Theta, multi = TRUE )
      output <- (append(output, dec)) 
    }
    return(output)
  } else if (method  == "EM"){
    output <- multi_HMM_EM( x = x, theta = theta, m = m, L1 = L1, L2 = L2, L3 = L3,
                        L4 = L4, L5 = L5, iterations = iterations,
                        DELTA = DELTA )
    
    #Decoding
    if (decoding==TRUE){
      dec <- decode( x = x, m = m, L1 = L1, L2 = L2, L3 = L3, L4 = L4, L5 = L5,
                   gamma = output$Gamma, delta = output$Delta,
                   theta = output$Theta, multi = TRUE )
      output <- (append(output, dec)) 
    }
    return(output)
  }else {
    warning("The supplied method is not available in the HMM function.
            Chooose between \"DM\" or \"EM\".", call. = TRUE)
  }
}
