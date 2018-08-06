#' Forward Probability function
#' 
#' @description This function calculates the forward probability of the HMM states. In the general 
#' literatur also revered to as alpha matrix
#' 
#' @param m number of Likelihoods
#' @param N length of the supplied dataset
#' @param sigma sigma vector  
#' @param Gamma Gamma matrix 
#' @param p vector of Likelihood probabilities of dataset
#' @param set index vector to align the p vector 
#' 
#' @return Returns the alpha as a matrix with the dimension m as columns 
#' and the datapoints t as rows 




alpha_function<-function( m, N, sigma, Gamma, p, set ){
  
  alpha<-matrix(, ncol=m, nrow = N)
  alpha[1,]<-t(sigma)%*%diag(c(p[set]))
  
  for (t in 2:N){
    alpha[t,] <- alpha[t-1,]%*%Gamma%*%diag(c(p[set+t-1]))
  }
  return(alpha)
}
