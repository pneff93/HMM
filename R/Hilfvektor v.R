#' v-function
#'
#' @description This function calculates the v-function which is part of the estimation step.
#' Each element in v aligns with the conditional expectation to reach state k from the previous
#' state j given your datapoint x(t). Due to the normalizing of the factors we also need to include
#' the corresponding weights into our function.
#'
#'
#' @param m number of likelihoods
#' @param N length of the supplied dataset
#' @param beta beta matrix
#' @param u output matrix of the u_function
#' @param Gamma Gamma matrix
#' @param p vector of likelihood probabilities of dataset
#' @param set index vector to align the p vector
#'

#' @return returns the coresponding values BUT due to further calculation steps
#' the m x m matrix of v-values for each datapoint t is
#' is returned as one row in the output matrix.
#'




v_function<-function( m, N, beta, p, c, Gamma, set,u){


  v<-matrix(,nrow=0, ncol=m*m)
  out<-matrix(, nrow=m, ncol=m)

  for (t in 2:N){
    pp<-diag(c(p[set+t-1]))

    for (i in 1:m){
      for(j in 1:m){
        out[i,j] <- (u[t-1,i]*Gamma[i,j]*pp[j,j]*beta[i,t])/(beta[i,t-1]*c[t-1])

      }
    }
    v<-c(v, c(t(out)))



  }
  v<-matrix(v, nrow = N-1, byrow = N)

  return(v)
}
