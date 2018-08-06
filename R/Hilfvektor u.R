#' u-function
#' 
#' @description This function calculates the u-function which is part of the estimation step. 
#' Each element in u aligns with the conditional expectation to 
#' reach state j given your datapoint x(t)
#' 
#' 
#' @param m number of Likelihoods
#' @param N length of the supplied dataset
#' @param alpha alpha matrix 
#' @param beta beta matrix 
#' @param L Likelihood of alpha*beta'

#' @return returns the coresponding u's
#' 





u_function<-function(m, N, alpha, beta, L){
  
  u<-matrix(, ncol=m, nrow = N)
  
  for (t in 1:N){
    for (j in 1:m){
      u[t,j]<-(alpha[t,j]%*%t(beta[j,t]))/L
    }
  }
  return(u)
}
