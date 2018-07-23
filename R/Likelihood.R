Likelihood<-function(alpha, beta){
  
  L<-alpha[1,]%*%beta[,1]
  return(L)
}
