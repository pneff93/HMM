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
