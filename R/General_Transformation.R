#' Transformation function - DM
#' 
#' @description Due to our single optimation problem we have to build constrains for Gamma and
#' delta to fullfill the requirements of probabilties. 

#' 
#' @param factor see details
#' @param m number of Likelihoods 
#' 
#' @details In the direct maximisation the nlme()- minimisation function can not directly implement
#' the constrains of the parameter values. This function ensures that the the estimated parameters 
#' of the direct optimisation still fullfill their requirements. These are that the probabilities are
#' between zero and one and that the rows of gamma (as well as delta) sum up to one.
#' For this transformation the logit model is used. 
#' 
#' Thus with the input factor containing the elements that determine delta and Gamma we need 
#' the following number of elements for each parameter:
#' 
#' delta vector (1 x m)  - (m-1)  elements required
#' Gamma matrix (m x m)  - m(m-1) elements required
#' Theta vector (1 x m)  - m      elements required
#' 
#' By this defintion the input factor vector has to contain the elements for delta, Gamma and Theta
#' in that order and has to have the dimension: (m+1)(m-1) + m
#' 
#'
#' 
#' @return returns a matrix with the delta,Gamma and theta matrix bound together (collumn wise)
#' 
#' 



trans <- function (factor,m){
  # Building the constrains: 
  #delta via logit transformation:
  #with delta[1] <- exp(factor[1])/(1+exp(factor[1]))
  
  delta <-c()
  div <- 1 + sum(exp(factor[1:(m-1)]))
  for (i in 1:(m-1)){
    delta[i] <- exp(factor[i])/div
  }
  delta[m] <- 1/div
  
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
  

  theta <- factor[((m-1)*(m+1)+1):length(factor)]
  return(list(delta,gamma,theta))}



