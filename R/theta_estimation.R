theta_estimation<-function( m, N, u, x, theta, L1, L2, L3=NULL, L4=NULL, L5=NULL ){
  
  l3<-0
  
  #For the first Likelihood
  l3 <- l3 + u[,1]%*%(log(L1(x,theta[1])))
  #For the second Likelihood
  l3 <- l3 +u[,2]%*%(log(L2(x,theta[2])))
  
  #If we have more than two Likelihoods we can use this algorithm 
  if (m>2){
    l3 <- l3 +u[,3]%*%(log(L3(x,theta[3]))) 
  } else if (m>3) {
    l3 <- l3 +u[,4]%*%(log(L4(x,theta[4])))
  } else if (m>4){
    l3 <- l3 +u[,5]%*%(log(L5(x,theta[5])))
  }
  return(l3*-1)
}

