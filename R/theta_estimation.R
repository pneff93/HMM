theta_estimation<-function( m, N, u, x, theta, L1, L2, L3=NULL, L4=NULL, L5=NULL ){
  
  l3<-0
  
  for (j in 1:m){
    if (j==1){
      pp<-L1
    } else if (j==2){
      pp<-L2
    } else if (j==3){
      pp<-L3
    } else if (j==4){
      pp<-L4
    } else if (j==5){
      pp<-L5
    }
    
    for (t in 1:N){
      l3<-l3+u[t,j]*log( pp(x[t], theta[j])) 
    }
  }
  return(l3*-1)
}

