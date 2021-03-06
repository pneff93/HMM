
################################################################################
############################## multi_HMM_DM function ###########################


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
#' @param theta initial parameters for the estimation of the likelihood 
#' parameters. See details for more information. 
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
#' This function is able to calculate with multiple Theta values for the 
#' individual likelihoods. For each likelihood the right number initial starting
#' parameter has to be set in order to compute the estimation of the 
#' corresponding Thetas.
#' For each Likelihood the starting values must be in the format of a vector, 
#' which is then saved as a list element. 
#' 
#' e.g.: theta[[i]] <- c(parameter1, parameter2, ...)
#' 
#' The function then extracts the right number of parameters per likelihood and 
#' optimizes the values. 
 


multi_HMM_DM <- function( x, theta,  m,  L1, L2, L3, L4, L5 ){
  
  #The definition of the likelihood is in a seperate R-file to reduce complexity
  #of the code. It can be found under multi_LH
  
  
  ##########################
  #Extraction of Thetas
  

  #Extract the right number of Thetas for each individual likelihood and save 
  #them in a combined vector theta_hat that will be maximized later

  #Additionaly we extract for each likelihood the length of the corresponding 
  #Theta vector and use it to build the index vector start_index. 
  #To later extract not only the starting of the next Theta but also the end of 
  #the corresponding Theta we add a dummy index as the last element of 
  #start_index. Thus e.g. to extract the start and end index of the second vector
  #we compute: 
  #start index: start_index[2]
  #end index: start_index[3]-1
  

  start_index <- c()
  theta_hat <- c()
  
  for (i in 1:length(theta)){
    start_index[i] <- length(theta_hat)+1
    theta_hat <- c(theta_hat, theta[[i]])
  }
  #dummy index 
  start_index[length(start_index)+1] <- length(theta_hat)+1
  
  
  ##########################
  #Setting starting values
  
  #Delta and Gamma
  #factor starting value of all Delta and Gamma values are 1/m.
  #Because we transform the starting values in the function LH we take 
  #the transformation into account with the reverse link function of the logit
  #model 
  
  factor <- c()
  factor[1:(m-1)] <- log((1/m) / (1-(m-1)*(1/m)))
  factor[m:((m+1)*(m-1))] <- (log((1/m) / (1-(m-1) * (1/m))))
  
  #Theta
  factor <- c(factor, theta_hat)
  
  
  ##########################
  #Maximisation
  
  #We maximize the log-likelihood with nlminb()
  #The multi_LH function can be found beneath this function
  #Depending on the number of likelihoods
  if (m==2){
    factor_out <- nlminb(start = factor, multi_LH, x = x, m = m, L1 = L1, L2 = L2,
                        start_index = start_index)$par
  } else if (m==3) {
    factor_out <- nlminb(start = factor, multi_LH, x = x, m = m, L1 = L1, L2 = L2,
                        L3 = L3, start_index = start_index)$par
  } else if (m==4) {
    factor_out<- nlminb(start = factor, multi_LH, x = x, m = m, L1 = L1, L2 = L2, 
                        L3 = L3, L4 = L4, start_index = start_index)$par
  } else if (m==5) {
    factor_out <- nlminb(start = factor, multi_LH, x = x, m = m, L1 = L1, L2 = L2,
                        L3 = L3, L4 = L4, L5 = L5,
                        start_index = start_index)$par
  }
  
  ##########################
  #Transformation and output
  
  #Transform the maximized values to our Gamma/Delta/Theta and return the output
  out <- trans(factor_out, m)
  theta_hat <- out[[3]]
  theta_out <- list()
  
  for (i in 1:m){
    theta_out[[i]] <- round(theta_hat[start_index[i]:(start_index[i+1]-1)], 3)
  }
  
  final <- list(
    "Method of estimation:" = "Direct Maximisation of the likelihoods",
    Delta = round( out[[1]], 3),
    Gamma = round(out[[2]], 3),
    Theta = theta_out)
  return(final)
}





################################################################################
############################## multi_LH function ###############################


#' Likelihood of the Hidden Markov Model
#'
#' @description This function calculates the log-likelihood of the HMM mode.
#' For this, the scaled forward probabilities are computed. In respect of the LH
#' function, the multi_LH function incorporates the multi Thetas for the 
#' indivudal likelihoods. 
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
#' 
#' @details This function computes the log-likelihood of the forward 
#' probabilities of the HMM. Given the fact that the inputed factor vector is 
#' not restricted, we need to apply a transformation to transform the factor
#' variables to our suitable canditates Delta, Gamma and Theta. For this we apply 
#' the function trans().
#' 
#' The likelihood is constructed using the scaling/normalizing of the forward 
#' probabilities such that for each time point the sum of the foward probability
#' is equal to one. The scaling is neccesary to tackle the underflow problem, 
#' that arrises with an increasing sample size.
#' 
#' The multi_LH function only differes in the computation of the likelihood 
#' probability vector p, because we imput multiple Thetas in the corresponding 
#' individual likelihoods. 



multi_LH <- function( factor, x, m, L1, L2, L3 = NULL, L4 = NULL, L5 = NULL,
                      start_index ){
  
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
  
  p1 <- L1(x, theta[start_index[1]:(start_index[2]-1)])
  p2 <- L2(x, theta[start_index[2]:(start_index[3]-1)])
  p <- c(p1, p2)
  
  if(!is.null(L3)){
    p3 <- L3(x, theta[start_index[3]:(start_index[4]-1)])
    p <- c(p, p3)
  }
  if(!is.null(L4)){
    p4 <- L4(x, theta[start_index[4]:(start_index[5]-1)])
    p <- c(p, p4)
  }
  if(!is.null(L5)){
    p5 <- L5(x, theta[start_index[5]:(start_index[6]-1)])
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