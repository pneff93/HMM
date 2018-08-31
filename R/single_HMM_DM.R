
################################################################################
############################## Single_HMM_DM function ##########################


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



HMM_DM<-function( x, m, L1, L2, L3, L4, L5 ){

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
  #The Log-Likelihood can be found below in the LH function 
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





################################################################################
############################## LH function #####################################


#' Likelihood of the Hidden Markov Model
#'
#' @description This function calculates the log-likelihood of the HMM mode.
#' For this, the scaled forward probabilities are computed. 
#'
#' @param factor Input of variables that are unrestricted
#' @param x  a sample of a Hidden Markov Model
#' @param m the number of states
#' @param L1 likelihood of the first hidden state
#' @param L2 likelihood of the second hidden state
#' @param L3 optional. likelihood of the third hidden state
#' @param L4 optional. likelihood of the 4th hidden state
#' @param L5 optional. likelihood of the 5th hidden state
#'
#' @return negative Likelihood
#' 
#' @details This function computes the log-likelihood of the forward 
#' probabilities of the HMM. Given the fact that the inputed factor vector is 
#' not restricted, we need to apply a transformation to transform the factor
#' variables to our suitable canditates Delta, Gamma and Theta. For this we apply 
#' the function trans().
#' 
#' The likelihood is constructed using the scaling/normalizing of the forward 
#' probabilities such that for each time point the sum of the foward probabily is
#' equal to one. The scaling is neccesary to tackle the underflow problem, that 
#' arrises with an increasing sample size.



LH <- function ( factor, x, m, L1, L2, L3 = NULL, L4 = NULL, L5 = NULL ){
  
  ##########################
  #Transformation 
  
  #We first have to transform the factors without constraints into our
  #Delta/Gamma/Theta values with constraints
  out <- trans(factor, m)
  delta <- out[[1]]
  gamma <- out[[2]]
  theta <- out[[3]]
  
  N <- length(x)
  
  #likelihoods
  
  #We combine the individual likelihoods for each timepoint to a general 
  #probability vector. The index set defines every cut, where a new likelihood
  #begins. This measure increases the performance of our calculations drastically
  #(instead of calculating a huge diagonal matrix)
  
  p1 <- L1(x, theta[1])
  p2 <- L2(x, theta[2])
  p <- c(p1,p2)
  
  if(!is.null(L3)){
    p3 <- L3(x, theta[3])
    p <- c(p, p3)
  }
  if(!is.null(L4)){
    p4 <- L4(x, theta[4])
    p <- c(p, p4)
  }
  if(!is.null(L5)){
    p5 <- L5(x, theta[5])
    p <- c(p, p5)
  }
  
  set <- seq(1, length(p)-N+1, length.out = m)
  
  
  
  ##########################
  #Computation of normalized forward probabilities
  
  #Computation of the log-likelihood with normalized alphas to tackle the 
  #underflow problem
  
  nalpha <- matrix(, nrow = m, ncol = N)
  v <- delta %*% diag(c(p[set]))
  u <- sum(v)
  l <- log(u)
  nalpha [, 1] <- t(v/u)
  
  
  for (t in 2:N){
    v <- nalpha[, t-1] %*% gamma %*% diag(c(p[set+t-1]))
    u <- sum(v)
    l <- l + log(u)
    nalpha[, t] <- t(v/u)
    
  }
  return (-1*l)
}