# Richtig harte Spielwiese 

#Transformation Function:

#' Fitting a Hidden Markov Model
#' @description Estimation of the transition probabilites, the initial state probabilites and the hidden state parameters of a Hidden Markov Model
#' by using the direct Maximisation of the Likelihoods
#' @param x, a sample of a Mixed Model
#' @param m, the number of states
#' @param L1, likelihood of the first hidden state
#' @param L2, likelihood of the second hidden state
#' @param L3-L5, optional. likelihoods of the third, 4th and 5th hidden state
#' @return The estimated parameters.
#' @details This function estimates the Hidden Markov states by maximising the normalized Log-Likelihood
#' of the forward propabilities. Due to the fact that bot the Gamma-matrix as well as the Sigma-matrix 
#' have some constrains, the function includes the restrictions within the function. 
#' 




HMM3<-function(x, m, L1, L2){
  
  #This function is the Log-Likelihood of the forward propabilities of the HMM
  #the input factor are the transformed values of Sigma,Gamma and Theta that hold
  # under the constrains
  
  LH <- function (factor,x,m, L1, L2){
    
    
    #We first have to transform the factors without constrains into our Sigma/Gamma/Theta 
    #values with constrains 
    out <- trans(factor,m)
    sigma <- out[,1]
    gamma <- out[,c(-1,-ncol(out))]
    theta <- out[,ncol(out)]
    
    T <- length(x)
    
    p1<-L1(x, theta[1])
    p2<-L2(x, theta[2])
    p<-c(p1,p2)
    set<-seq(1, length(p)-T+1, length.out = m)
    
    #Computation of the Log-Likelihood with normalized alphas to tackle the underflow problem 
 
    #normalized alpha
    nalpha<- matrix(,nrow=m, ncol=T)
    v <- sigma%*%diag(c(p[set]))
    u <- sum(v)
    l <- log(u)
    nalpha [,1] <- t(v/u)
    
 
    for (t in 2:T){
      v <- nalpha[,t-1]%*%gamma%*%diag(c(p[set+t-1]))
      u <- sum(v)
      l <- l + log(u)
      nalpha[,t]<- t(v/u)
      
    }
    return (-1*l)
  }
  
  #setting starting values: 
  
  #factor starting value (sigma and all gamma values are 1/m)
  #e.g for m=2 sigma= (0.5,0.5) 
  #and Gamma =(0.5,0.5)
  #           (0.5,0.5)
  #due to the transformation we have to use the reverse link function of the probit model 
  factor <- c()
  factor[1:(m-1)] <- log((1/m)/(1-(m-1)*(1/m)))
  factor[m:((m+1)*(m-1))] <-(log((1/m)/(1-(m-1)*(1/m))))
  factor[((m-1)*(m+1)+1):((m-1)*(m+1)+m)]<-sample(x, size=m)
  factor
  
  
  
  
  #Now we maximize the log-Likelihood with the nlminb 
   factor_out <- nlminb(start=factor,LH,x=x,m=m,L1=L1,L2=L2)$par
   
   
  #Transform the maximized values to our Gamma/Sigma/Theta and return the output 
   out <- trans(factor_out,m)
   final <- list("Sigma Value" = sigma <- out[,1],
        "Gamma Matrix" = gamma <- out[,c(-1,-ncol(out))],
        "Estimated Theta of the m-Likelihoods" = theta <- out[,ncol(out)])
   
   return(final)
}  
  