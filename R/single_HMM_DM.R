#' Fitting a Hidden Markov Model via the Direct Maximisation
#'
#' @description Estimation of the transition probabilites,  
#' the initial state probabilites and the hidden state parameters of a Hidden 
#' Markov Model by using the Direct Maximisation of the global log-likelihood.
#'
#'
#' @param x a sample of a Hidden Markov Model
#' @param m the number of states
#' @param L1 likelihood of the first hidden state
#' @param L2 likelihood of the second hidden state
#' @param L3 optional. likelihood of the third hidden state
#' @param L4 optional. likelihood of the 4th hidden state
#' @param L5 optional. likelihood of the 5th hidden state
#'
#' @return The estimated parameters are rounded by 3 decimals and returned in a 
#' list.
#' 
#' @details This function estimates the Hidden Markov states by maximising the 
#' normalized log-likelihood of the forward propabilities. Due to the fact that
#' both the Gamma matrix as well as the Delta vector have some constraints, the 
#' function first applies some restrictions and then uses the base-R maximisation
#' to gain the most likely variables. 
#'


HMM3<-function( x, m, L1, L2, L3, L4, L5 ){

  #The definition of the likelihood is in a seperate R-file to reduce complexity
  #of the code. It can be found under single_LH.R
  
  
  ##########################
  #Setting starting values:
  
  #Delta and Gamma 
  #factor starting value of all Delta and Gamma values are 1/m.
  #Because we transform the starting values in the function LH we  take the 
  #the transformation into account with the reverse link function of the logit
  #model 

  factor <- c()
  factor[1:(m-1)] <- log((1/m) / (1-(m-1)*(1/m)))
  factor[m:((m+1)*(m-1))] <- (log((1/m) / (1-(m-1)*(1/m))))
  
  #Theta
  #The starting values of the Theta values are the quantiles of the distribution
  #e.g. for m = 2 we use q_1 = 0.333 and q_2 = 0.666 as quantile values.

  #quantile defintion vector
  s <- seq(0, 1, length.out = m+2)[c(-1, -(m+2))]
  factor[((m-1)*(m+1)+1):((m-1)*(m+1)+m)] <- quantile(x, s)

  
  ##########################
  #Maximisation 
  
  #We maximize the log-likelihood with nlminb()
  #Depending on the number of likelihoods
  if (m==2){
   factor_out <- nlminb(start = factor, LH, x = x, m = m, L1 = L1, L2 = L2)$par
  } else if (m==3) {
  factor_out <- nlminb(start = factor, LH, x = x, m = m, L1 = L1, L2 = L2,
                      L3 = L3)$par
  } else if (m==4) {
  factor_out <- nlminb(start = factor, LH, x = x, m = m, L1 = L1, L2 = L2,
                      L3 = L3, L4 = L4)$par
  } else if (m==5) {
  factor_out<- nlminb(start = factor, LH, x = x, m = m, L1 = L1, L2 = L2, 
                      L3 = L3, L4 = L4, L5 = L5)$par
  }

  
  ##########################
  #Transformation and output
  
  #Transform the maximized values to our Gamma/Delta/Theta and return the output
   out <- trans(factor_out, m)
   
   final <- list(
        "Method of estimation:" = "Direct Maximisation of the likelihoods", 
        Delta = round( out[[1]], 3), 
        Gamma = round(out[[2]], 3), 
        Theta = round( out[[3]], 3))
   return(final)
}
