beta_function<-function( m, N, Gamma, p, set ){
  
  beta<-matrix(, ncol = N, nrow = m)
  beta[,N]<-rep(1, times=m)
  
  for (t in (N-1):1){
    beta[,t] <- Gamma%*%diag(c(p[set+t]))%*%beta[,t+1]
  }
  
  return(beta)
}
