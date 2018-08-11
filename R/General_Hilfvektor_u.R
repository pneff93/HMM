#' Normalized u-function
#'
#' @description This function calculates the u-function which is part of the estimation step.
#' Each element in u aligns with the conditional expectation to
#' reach state j given your datapoint x(t). The corresponding weights are already included in the
#' alpha and beta values.
#'
#'
#' @param m number of likelihoods
#' @param N length of the supplied dataset
#' @param alpha normalized alpha matrix
#' @param beta normalized beta matrix
#'

#' @return returns the corresponding normalized u's
#'


u_function<-function(m, N, alpha, beta){


  u<-matrix(, ncol=m, nrow = N)

  for (t in 1:N){
    for (j in 1:m){
      u[t,j]<-(alpha[t,j]%*%t(beta[j,t]))/(alpha[t,]%*%beta[,t])
    }
  }
  return(u)
}

