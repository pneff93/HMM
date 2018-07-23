u_function<-function(m, N, alpha, beta, L){
  
  u<-matrix(, ncol=m, nrow = N)
  
  for (t in 1:N){
    for (j in 1:m){
      u[t,j]<-(alpha[t,j]%*%t(beta[j,t]))/L
    }
  }
  return(u)
}
