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
#'
#' @return The estimated parameters are rounded by 3 decimals and returned in a list.
#' @details This function estimates the Hidden Markov states by maximising the normalized log-likelihood
#' of the forward propabilities. Due to the fact that bot the Gamma-matrix as well as the Sigma-matrix
#' have some constraints, the function includes the restrictions within the function.
#'


HMM3<-function(x, m, L1, L2,L3,L4,L5){

  #This function is the log-likelihood of the forward propabilities of the HMM
  #the input factor are the transformed values of Sigma, Gamma and Theta that hold
  #under the constrains

  LH <- function (factor,x,m, L1, L2,L3=NULL,L4=NULL,L5=NULL){


    #We first have to transform the factors without constrains into our Sigma/Gamma/Theta
    #values with constrains
    out <- trans(factor,m)
    sigma <- out[,1]
    gamma <- out[,c(-1,-ncol(out))]
    theta <- out[,ncol(out)]

    T <- length(x)


    #likelihoods
    p1<-L1(x, theta[1])
    p2<-L2(x, theta[2])
    p<-c(p1,p2)

    if(!is.null(L3)){
      p3<-L3(x, theta[3])
      p<-c(p,p3)
    }
    if(!is.null(L4)){
      p4<-L4(x, theta[4])
      p<-c(p,p4)
    }
    if(!is.null(L5)){
      p5<-L5(x, theta[5])
      p<-c(p,p5)
    }

    set<-seq(1, length(p)-T+1, length.out = m)

    #Computation of the log-likelihood with normalized alphas to tackle the underflow problem

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
  #due to the transformation we have to use the reverse link function of the probit model

  factor <- c()
  factor[1:(m-1)] <- log((1/m)/(1-(m-1)*(1/m)))
  factor[m:((m+1)*(m-1))] <-(log((1/m)/(1-(m-1)*(1/m))))
  
  theta_hat<-c()
  #quantile defintion vector
  s <- seq(0,1,length.out = m+2)[c(-1,-(m+2))]
  
  factor[((m-1)*(m+1)+1):((m-1)*(m+1)+m)]<-quantile(x,s)



  #We maximize the log-likelihood with nlminb
  #Depending on the number of likelihoods
  if (m==2){
   factor_out<- nlminb(start=factor,LH,x=x,m=m,L1=L1,L2=L2)$par
  } else if (m==3) {
  factor_out<- nlminb(start=factor,LH,x=x,m=m,L1=L1,L2=L2,L3=L3)$par
  } else if (m==4) {
  factor_out<- nlminb(start=factor,LH,x=x,m=m,L1=L1,L2=L2,L3=L3,L4=L4)$par
  } else if (m==5) {
  factor_out<- nlminb(start=factor,LH,x=x,m=m,L1=L1,L2=L2,L3=L3,L4=L4,L5=L5)$par
  }


  #Transform the maximized values to our Gamma/Sigma/Theta and return the output
   out <- trans(factor_out,m)
   final <- list(
        "Method of estimation:" = "Direct Maximisation of the likelihoods",
        Sigma = round( out[,1],3),
        Gamma= round(out[,c(-1,-ncol(out))],3),
        Theta = round( out[,ncol(out)],3))
   return(final)
}
