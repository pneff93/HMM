alpha_function<-function( m, N, sigma, Gamma, p, set ){

  alpha<-matrix(, ncol=m, nrow = N)
  alpha[1,]<-t(sigma)%*%diag(c(p[set]))

  for (t in 2:N){
    prod<-diag(1, m)
    for (s in 2:t){
      prod<-prod%*%Gamma%*%diag(c(p[set+s-1]))
    }
    alpha[t,]<-alpha[1,]%*%prod
  }
  return(alpha)
}
