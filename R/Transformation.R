#Due to our single optimation problem we have to build constrains for Gamma, Sigma and 
# maybe even lambda
#to tackle this problem we let the factors unconstrained and implement the constrains within
#the function. After optimization we then rebuild the constrains

#############
#Der Factor ist ein Vector, die alle drei variablen collumwise gestapelt besitzt.
#also cbind(sigma,gamma,theta), wobei theta die einzelenen Faktoren der Likelihoods 
#als Vektor sind. die Dimensionen 1x (m+m*m+m)

trans <- function (factor,m){
  # Building the constrains: 
  #sigma via probit transformation:
  #with sigma[1] <- exp(factor[1])/(1+exp(factor[1]))
  
  sigma <-c()
  div <- 1 + sum(exp(factor[1:(m-1)]))
  for (i in 1:(m-1)){
    sigma[i] <- exp(factor[i])/div
  }
  sigma[m] <- 1/div
  
  #Building restrictions on Gamma
  #for i =! j.
  #As example for m=3 
  # ( --- tau12 tau13)
  # (tau21 ---  tau23)
  # (tau31 tau32 ---)
  
  #we extract therefore the next m(m-1) elements of factor 
  #so the elements:  (m-1)+1=m till (m-1)+m(m-1)=(m+1)(m-1)
  x <- matrix(factor[m:((m+1)*(m-1))],ncol=(m-1),nrow=m,byrow=T)
  
  
  gamma <- matrix(,nrow=m,ncol=m)
  
  for (i in 1:m){
    div <- 1+sum(exp(x[i,]))
    count <- 0
    
    for (j in 1:m){
      if(i==j){ 
        gamma[i,j] <- 1/div
      } else {
        count <- count + 1
        #we need a counter to asign the right values from x 
        gamma[i,j]<- exp(x[i,count])/div
      }
    }
  }
  theta <- factor[((m-1)*(m+1)+1):((m-1)*(m+1)+m)]
  return(cbind(sigma,gamma,theta))}



