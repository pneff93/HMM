alpha_function<-function( m, N, sigma, Gamma, p, set ){
  
  alpha<-matrix(, ncol=m, nrow = N)
  alpha[1,]<-t(sigma)%*%diag(c(p[set]))
  
  for (t in 2:N){
    alpha[t,] <- alpha[t-1,]%*%Gamma%*%diag(c(p[set+t-1]))
  }
  return(alpha)
}
