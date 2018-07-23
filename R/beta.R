beta_function<-function( m, N, Gamma, p, set ){

  beta<-matrix(, ncol = N, nrow = m)
  beta[,N]<-rep(1, times=m)

  for (t in 1:(N-1)){
    prod<-diag(1, m)
    for (s in (t+1):N){
      prod<-prod%*%Gamma%*%diag(c(p[set+s-1]))
    }
    beta[,t]<-prod%*%rep(1, times=m)
  }
  return(beta)
}
