#' Fitting a Hidden Markov Model via the Direct Maximisation
#'
#' @description Estimation of the transition probabilites, the initial state probabilites and the hidden state parameters of a Hidden Markov Model
#' by using the Direct Maximisation of the likelihoods.
#'
#'
#' @param x a sample of a Mixed Model
#' @param m the number of states
#' @param L1 likelihood of the first hidden state
#' @param L2 likelihood of the second hidden state
#' @param L3 optional. likelihood of the third hidden state
#' @param L4 optional. likelihood of the 4th hidden state
#' @param L5 optional. likelihood of the 5th hidden state
#' @param theta initial parameters for the estimation of the likelihood parameters 
#' @return The estimated parameters are rounded by 3 decimals and returned in a list.
#' @details This function estimates the Hidden Markov states by maximising the normalized log-likelihood
#' of the forward propabilities. Due to the fact that bot the Gamma-matrix as well as the Sigma-matrix
#' have some constraints, the function includes the restrictions within the function.
#'


multiHMM3<-function(x,theta, m, L1, L2,L3,L4,L5){
  
  #The calulation of the multi LH function is in a seperate function call 
  #setting starting values:
  
  #factor starting value (sigma and all gamma values are 1/m)
  #due to the transformation we have to use the reverse link function of the probit model
  
  #Like in the EM-function we first decode the theta list into a theta vector
  #this theta hat is later crutial as initial parameter for the optimisation
  #set2 is the index that identifies the next set of theta parameters
  set2 <- c()
  theta_hat <-c()
  
  for (i in 1:length(theta)){
    set2[i] <- length(theta_hat)+1
    theta_hat <-c(theta_hat,theta[[i]])
  }
  #We add one more index to know, where the last one would end
  #(or the additional start)
  set2[length(set2)+1] <- length(theta_hat)+1
  
  
  #computation of the initial factors 
  
  factor <- c()
  factor[1:(m-1)] <- log((1/m)/(1-(m-1)*(1/m)))
  factor[m:((m+1)*(m-1))] <-(log((1/m)/(1-(m-1)*(1/m))))
  
  factor <- c(factor,theta_hat)
  
  
  
  #We maximize the log-likelihood with nlminb
  #Depending on the number of likelihoods
  if (m==2){
    factor_out<- nlminb(start=factor,multiLH,x=x,m=m,L1=L1,L2=L2,set2=set2)$par
  } else if (m==3) {
    factor_out<- nlminb(start=factor,multiLH,x=x,m=m,L1=L1,L2=L2,L3=L3,set2=set2)$par
  } else if (m==4) {
    factor_out<- nlminb(start=factor,multiLH,x=x,m=m,L1=L1,L2=L2,L3=L3,L4=L4,set2=set2)$par
  } else if (m==5) {
    factor_out<- nlminb(start=factor,multiLH,x=x,m=m,L1=L1,L2=L2,L3=L3,L4=L4,L5=L5,set2=set2)$par
  }
  
  
  #Transform the maximized values to our Gamma/Sigma/Theta and return the output
  out <- trans(factor_out,m)
  theta_hat <- out[[3]]
  theta_out <-list()
  
  for (i in 1:m){
    theta_out[[i]] <- round(theta_hat[set2[i]:(set2[i+1]-1)],3)
  }
  
  final <- list(
    "Method of estimation:" = "Direct Maximisation of the likelihoods",
    Sigma = round( out[[1]],3),
    Gamma= round(out[[2]],3),
    Theta = theta_out)
  return(final)
}
