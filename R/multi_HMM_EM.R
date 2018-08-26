#' Fitting a Hidden Markov Model via Expectation-Maximisation Algorithm
#'
#' @description Estimation of the transition probabilites, the initial state probabilites and the hidden state parameters of a Hidden Markov Model
#' by using the Baum-Welch Algorithm.
#'
#' @param x a sample of a Mixed Model
#' @param m the number of states
#' @param L1 likelihood of the first hidden state
#' @param L2 likelihood of the second hidden state
#' @param L3 optional. likelihood of the third hidden state
#' @param L4 optional. likelihood of the 4th hidden state
#' @param L5 optional. likelihood of the 5th hidden state
#' @param theta initial parameters for the estimation of the likelihood parameters. See details for more information 
#' @param iterations optional. number of iterations for the EM-Algorithm
#' @param DELTA optional. stop criterion for the EM-Algorithm
#'
#' @return The estimated parameters are rounded by 3 decimals and returned in a list.
#' @details This function estimates the hidden states of the Hidden Markov Model with the
#' help of the Baum Welch algorithm. The function iteratively applies both estimation and
#' maximisation steps to arrive at the predicted parameters. When the maximum difference between
#' the functions is abitrarily small (see DELTA) the iteration stops.
#'
#' For each individual likelihood an initial starting parameter has to be set in order to compute the estimation of the corresponding
#' Thetas. Each groupe of parameters is placed in a seperate element of the theta list as a vector e.g.:
#' theta[[i]] <- c(parameter1,parameter2,...)



multiHMM2<-function(x,theta, m, L1, L2, L3=NULL, L4=NULL, L5=NULL, iterations=NULL, DELTA=NULL){
  
  #setting starting values:
  
  theta_hat<-c()
  set2 <- c()
  
  #Here the theta list gets transformed into a vector for 
  #further maximisation
  #The vector set2 is also added which is an index, where the
  #new parameters start.
  for (i in 1:length(theta)){
    set2[i] <- length(theta_hat)+1
    theta_hat <-c(theta_hat,theta[[i]])
  }
  #We add one more index to know, where the last one would end
  #(or the additional start)
  set2[length(set2)+1] <- length(theta_hat)+1
  
  Gamma_hat<-matrix(, ncol = m, nrow = m)
  Gamma_hat[,]<-1/m
  
  delta_hat<-matrix(, nrow = m)
  delta_hat[,]<-1/m
  
  
  #sample size
  N<-length(x)
  
  #likelihoods
  
  
  
  #We define a treshold of the sample size for every likelihood
  set<-seq(1, N*m-N+1, length.out = m)
  
  
  #while-loop parameters
  if(is.null(iterations)){
    z<-500
  } else {
    z<-iterations
  }
  
  if(is.null(DELTA)){
    d<-0.0001
  } else {
    d<-DELTA
  }
  
  q<-1
  DELTA<-Inf
  
  
  while (DELTA>d && q<z){
    
    q<-q+1
    
    theta<-theta_hat
    Gamma<-Gamma_hat
    delta<-delta_hat
    
    #computation of the individual likelihoods of each state 
    
    p1<-L1(x, theta_hat[set2[1]:(set2[2]-1)])
    p2<-L2(x, theta_hat[set2[2]:(set2[3]-1)])
    p<-c(p1,p2)
    
    if(!is.null(L3)){
      p3<-L3(x, theta_hat[set2[3]:(set2[4]-1)])
      p<-c(p,p3)
    }
    if(!is.null(L4)){
      p4<-L4(x, theta_hat[set2[4]:(set2[5]-1)])
      p<-c(p,p4)
    }
    if(!is.null(L5)){
      p5<-L5(x, theta_hat[set2[5]:(set2[6]-1)])
      p<-c(p,p5)
    }    
   
    
    
    out<-alpha_function(m, N, delta, Gamma, p, set)
    alpha <- out[[1]]
    weight <- out[[2]]
    
    beta<-beta_function(m, N, Gamma, p, set, c= weight )
    
    
    u<-u_function(m, N, alpha, beta)
    
    v<-v_function(m, N,beta, p,c=weight, Gamma, set, u = u)
    
    
    #set in the likelihood (consisting of three parts)
    delta_hat<-delta_hat<-u[1,]
    
    f<-apply(v,2,sum)
    f<-matrix(f, ncol = m, byrow = N)
    
    Gamma_hat<-matrix(, ncol = m, nrow = m)
    
    for (j in 1:m){
      Gamma_hat[j,]<-f[j,]/sum(f[j,])
    }
    
    #We maximize the third part of the likelihood with nlminb
    #Depending on the number of likelihoods
    
    if(m==2){
      theta_hat<-nlminb(start = theta, multi_theta_estimation, m=m, N=N, u=u, x=x, L1=L1, L2=L2, set2=set2)$par
    } else if (m==3){
      theta_hat<-nlminb(start = theta, multi_theta_estimation, m=m, N=N, u=u, x=x, L1=L1, L2=L2, L3=L3, set2=set2)$par
    } else if (m==4){
      theta_hat<-nlminb(start = theta, multi_theta_estimation, m=m, N=N, u=u, x=x, L1=L1, L2=L2, L3=L3, L4=L4, set2=set2)$par
    } else if (m==5){
      theta_hat<-nlminb(start = theta, multi_theta_estimation, m=m, N=N, u=u, x=x, L1=L1, L2=L2, L3=L3, L4=L4, L5=L5, set2=set2)$par
      
    }
    
    #DELTA Estimation (do not confuse with delta-vector)
    DELTA1<-max(abs(delta_hat-delta))
    DELTA2<-max(abs(Gamma_hat-Gamma))
    DELTA3<-max(abs(theta_hat-theta))
    DELTA<-max(DELTA1, DELTA2, DELTA3)
    
  }
  
  
  #Build the output 
  theta_out <- list()
  for (i in 1:m){
    theta_out[[i]] <- round(theta_hat[set2[i]:(set2[i+1]-1)],3)
  }
  
  return<-list( "Method of estimation:" = "Maximisation via EM-Algorithm",
                delta=round(delta_hat,3),
                Gamma=round(Gamma_hat,3),
                Theta=theta_out,
                Iterations=q,
                DELTA=DELTA )
  
  return(return)
}

