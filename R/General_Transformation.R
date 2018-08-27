#' Transformation function - DM
#' 
#' @description The transformation function transforms a non-restricted parameter
#' vector into a restricted Delta, Gamma and Theta output.  
#' 
#' @param factor vector of unrestricted parameters
#' @param m number of Likelihoods 
#' 
#' @details In the direct maximisation the nlminb()- minimisation function can not 
#' directly implement the constraints of the parameter values. 
#' This function ensures that the estimated parameters of the direct 
#' optimisation still fullfill their requirements. The requirements are that all
#' probabilities are between zero and one and that the rows of Gamma (as well as
#' the vector Delta) sum up to one. For this transformation the logit model is 
#' used. 
#' 
#' From the input vector factor, the elements are extracted as follows and used 
#' for transformation: 
#' 
#' Delta vector (1 x m): The logit transformation requires (m-1)  elements, thus 
#' the first (m-1) elements of the factor vector are used. 
#'  
#' Gamma matrix (m x m): The logit transformation requires per row (m-1) elements
#' thus in total m(m-1) elements required. Accordingly the elements from index 
#' m til (m+1)(m-1) are extracted from the factor vector. 
#' 
#' Theta vector: Due to the fact that the trans() function can be used both from 
#' the single and multifactor Direct Maximisation the Theta vector has no pre-
#' defined lenght, but needs to have at least the length of m (thus each 
#' likelihood has at least one parameter). 
#' 
#' The resulting Delta, Gamma and Theta values are returned as a single list. 
#' 
#' 
#' 
#' 
#' @return List of Delta, Gamma, Theta
#' 
#' 
#' 



trans <- function( factor, m ){
  ##########################
  #Transformation Delta-vector
  
  #Building the constrains for Delta via logit transformation:
  # delta[i] <- exp(factor[i])/sum(1+exp(factor)) for i = 1 - (m-1)
  # delta[m] <- 1/sum(1+exp(factor))
  
  delta <- c()
  div <- 1 + sum(exp(factor[1:(m-1)]))
  for (i in 1:(m-1)){
    delta[i] <- exp(factor[i]) / div
  }
  delta[m] <- 1 / div
  
  
  ##########################
  #Transformation Gamma-matrix
  
  #Building restrictions on Gamma we extract therefore the next m(m-1) elements
  #of factor; so the elements: (m-1)+1 = m till (m-1)+m(m-1) = (m+1)(m-1)
 
  #For the Gamma matrix we determine the off-diagonal elements: 
  
  #for i =! j: 
  #exp(factor[i])/sum(1+exp(factor))
  #for i == j 
  #1/sum(1+exp(factor))
  
  #As example for m = 3 
  # ( --- tau12 tau13)
  # (tau21 ---  tau23)
  # (tau31 tau32 ---)

  
  #Recall that each rowsum has to be equal to one 
  #The variable count ensures that only the off-diagonal elements are calculated
  #with the factor vector 
  
  x <- matrix(factor[m:((m+1)*(m-1))], ncol = (m-1), nrow = m, byrow = TRUE)
  
  
  gamma <- matrix(, nrow = m, ncol = m)
  
  for (i in 1:m){
    div <- 1+sum(exp(x[i, ]))
    count <- 0
    
    for (j in 1:m){
      if(i==j){ 
        gamma[i, j] <- 1/div
      } else {
        count <- count + 1
        gamma[i, j] <- exp(x[i, count])/div
      }
    }
  }
  
  ##########################
  #Transformation Theta-values
  
  
  #The Theta values are not manipulated and are just stored in an individual 
  #vector 
  
  theta <- factor[((m-1)*(m+1)+1):length(factor)]
  
  
  return(list(delta, gamma, theta))}



