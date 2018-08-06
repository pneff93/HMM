#' v-function
#' 
#' @description This function calculates the v-function which is part of the estimation step. 
#' Each element in v aligns with the conditional expectation to 
#' reach state k from the previous state j given your datapoint x(t)
#' 
#' 
#' @param m number of Likelihoods
#' @param N length of the supplied dataset
#' @param alpha alpha matrix 
#' @param beta beta matrix 
#' @param L Likelihood of alpha*beta'
#' @param Gamma Gamma matrix 
#' @param p vector of Likelihood probabilities of dataset
#' @param set index vector to align the p vector 
#' 

#' @return returns the coresponding values BUT due to further calculation steps
#' the m x m matrix of v-values for each datapoint t is 
#' is returned as one row in the output matrix. 
#' 




v_function<-function( m, N, alpha, beta, p, L, Gamma, set){

  v<-matrix(,nrow=0, ncol=m*m)
  out<-matrix(, nrow=m, ncol=m)

  for (t in 2:N){
    pp<-diag(c(p[set+t-1]))

    for (j in 1:m){
      for (k in 1:m){
        out[j,k]<-alpha[t-1,j]*Gamma[j,k]*pp[k,k]*beta[k,t]/L
      }
    }
    v<-c(v, c(t(out)))
  }
  v<-matrix(v, nrow = N-1, byrow = N)
  return(v)
}
