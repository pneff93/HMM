
################################################################################
############################## multi_HMM_EM function ###########################


#' Fitting a Hidden Markov Model via Expectation-Maximisation Algorithm
#'
#' @description Estimation of the transition probabilites, the initial state
#' probabilites and the hidden state parameters of a Hidden Markov Model
#' by using the Baum-Welch Algorithm.
#'
#' @param x a sample of a Hidden Markov Model
#' @param m the number of states
#' @param L1 likelihood of the first hidden state
#' @param L2 likelihood of the second hidden state
#' @param L3 optional. likelihood of the third hidden state
#' @param L4 optional. likelihood of the 4th hidden state
#' @param L5 optional. likelihood of the 5th hidden state
#' @param theta initial parameters for the estimation of the likelihood 
#' parameters. See details for more information 
#' @param iterations optional. number of iterations for the EM-Algorithm
#' @param DELTA optional. stop criterion for the EM-Algorithm
#'
#' @return The estimated parameters are rounded by 3 decimals and returned in a 
#' list.
#' 
#' @details This function estimates the hidden states of the Hidden Markov Model
#' with the help of the Baum Welch algorithm. The function iteratively applies 
#' both estimation and maximisation steps to arrive at the predicted parameters.
#' When the maximal difference between present and prior parameter is abitrarily 
#' small (defined with DELTA) the iteration stops.
#'
#' This function is able to calculate with multiple Theta values for the 
#' individual likelihoods. For each likelihood the right number initial starting
#' parameter has to be set in order to compute the estimation of the 
#' corresponding Thetas.
#' For each Likelihood the starting values must be in the format of a vector, 
#' which is then saved as a list element. 
#' 
#' e.g.: theta[[i]] <- c(parameter1, parameter2, ...)
#' 
#' The function then extracts the right number of parameters per likelihood and 
#' optimizes the values. 



multi_HMM_EM <- function( x, theta, m, L1, L2, L3 = NULL, L4 = NULL, L5 = NULL,
                    iterations = NULL, DELTA = NULL ){
  

  
  ##########################
  #Initiating starting values 
  
  #Extraction of Thetas
  #Extract the right number of Thetas for each individual likelihood and save 
  #them in a combined vector theta_hat, that will be maximized later
  
  #Additionaly we extract for each likelihood the length of the corresponding 
  #Theta vector and use it to build the index vector start_index. 
  #To later extract not only the starting of the next Theta but also the end of 
  #the corresponding Theta we add a dummy index as the last element of 
  #start_index. Thus e.g. to extract the start and end index of the second vector
  #we compute: 
  #start index: start_index[2]
  #end index: start_index[3]-1
  
  
  theta_hat <- c()
  start_index <- c()
  
  for (i in 1:length(theta)){
    start_index[i] <- length(theta_hat)+1
    theta_hat <- c(theta_hat, theta[[i]])
  }
  
  #dummy index
  start_index[length(start_index)+1] <- length(theta_hat)+1
  
  
  #Gamma
  Gamma_hat <- matrix(, ncol = m, nrow = m)
  Gamma_hat[, ]<-1/m
  
  #Delta
  delta_hat<-matrix(, nrow = m)
  delta_hat[, ]<-1/m
  
  
  #sample size
  N<-length(x)
  

  
  #We define a treshold of the sample size for every likelihood
  #In the iteration we combine the individual likelihoods to a general 
  #probability vector.
  set <- seq(1, N*m-N+1, length.out = m)
  
  ##########################
  #while-loop parameters
  if(is.null(iterations)){
    z <- 500
  } else {
    z <- iterations
  }
  
  if(is.null(DELTA)){
    d <- 0.0001
  } else {
    d <- DELTA
  }
  
  q <- 1
  DELTA <- Inf
  
  
  ##########################
  ##########################
  #Iterations 
  
  while (DELTA>d && q<z){
    
    q <- q+1
    
    theta <- theta_hat
    Gamma <- Gamma_hat
    delta <- delta_hat
    
    #computation of the individual likelihoods of each state 
    
    #We combine the individual likelihoods for each timepoint to a general 
    #probability vector. The index set defines every cut, where a new likelihood
    #begins. This measure increases the performance of our calculations 
    #drastically
    
    #Recall that the start_index can show a start and beginning of the
    #corresponding Theta value 
    
    p1 <- L1(x, theta_hat[start_index[1]:(start_index[2]-1)])
    p2 <- L2(x, theta_hat[start_index[2]:(start_index[3]-1)])
    p <- c(p1, p2)
    
    if(!is.null(L3)){
      p3 <- L3(x, theta_hat[start_index[3]:(start_index[4]-1)])
      p <- c(p, p3)
    }
    if(!is.null(L4)){
      p4 <- L4(x, theta_hat[start_index[4]:(start_index[5]-1)])
      p <- c(p, p4)
    }
    if(!is.null(L5)){
      p5 <- L5(x, theta_hat[start_index[5]:(start_index[6]-1)])
      p <- c(p, p5)
    }    
   
    ##########################
    #Calculation of forward and backward probability 
    
    out <- alpha_function( m, N, delta, Gamma, p, set )
    alpha <- out[[1]]
    weight <- out[[2]]
    
    beta <- beta_function( m, N, Gamma, p, set, weight = weight )
    
    
    ##########################
    #Estimation of u and v 
    
    u <- u_function( m, N, alpha, beta )
    
    v <- v_function( m, N, beta, p, weight = weight, Gamma, set, u = u )
    
    
    
    ##########################
    #Maximisation
    
    #Delta
    delta_hat <- delta_hat <- u[1, ]
    
    #Gamma
    f <- apply(v, 2, sum)
    f <- matrix(f, ncol = m, byrow = N)
    
    Gamma_hat <- matrix(, ncol = m, nrow = m)
    
    for (j in 1:m){
      Gamma_hat[j, ] <- f[j, ] / sum(f[j, ])
    }
    
    #Theta
    #We maximize the third part of the likelihood with nlminb()
    #Depending on the number of likelihoods
    #The multi_theta_estimation function is below this function in this file 
    
    if(m==2){
      theta_hat <- nlminb(start = theta, multi_theta_estimation, m = m, N = N,
                        u = u, x = x, L1 = L1, L2 = L2,
                        start_index = start_index)$par
    } else if (m==3){
      theta_hat <- nlminb(start = theta, multi_theta_estimation, m = m, N = N,
                        u = u, x = x, L1 = L1, L2 = L2, L3 = L3,
                        start_index = start_index)$par
    } else if (m==4){
      theta_hat <- nlminb(start = theta, multi_theta_estimation, m = m, N = N,
                        u = u, x = x, L1 = L1, L2 = L2, L3 = L3,L4 = L4,
                        start_index = start_index)$par
    } else if (m==5){
      theta_hat <- nlminb(start = theta, multi_theta_estimation, m = m, N = N,
                        u = u, x = x, L1 = L1, L2 = L2, L3 = L3, L4 = L4,
                        L5 = L5, start_index = start_index)$par
      
    }
    
    ##########################
    #DELTA Estimation (do not confuse with delta-vector)
    
    #We compute the difference between the estimated parameter and the prior 
    #parameter as a stopping condition for the loop. 
    DELTA1 <- max(abs(delta_hat-delta))
    DELTA2 <- max(abs(Gamma_hat-Gamma))
    DELTA3 <- max(abs(theta_hat-theta))
    DELTA <- max(DELTA1, DELTA2, DELTA3)
    
  }
  
  ##########################
  #Build the output 
  theta_out <- list()
  for (i in 1:m){
    theta_out[[i]] <- round(theta_hat[start_index[i]:(start_index[i+1]-1)], 3)
  }
  
  return <- list( "Method of estimation:" = "Maximisation via EM-Algorithm",
                Delta = round(delta_hat, 3),
                Gamma = round(Gamma_hat, 3),
                Theta = theta_out,
                Iterations = q,
                DELTA = DELTA )
  
  return(return)
}





################################################################################
############################## multi_theta_estimation_function #################


#' Multi Theta Estimation function
#'
#' @description This function calculates part of the global log-likelihood that 
#' is only dependent on the Theta value. Due to its proportionality, it is 
#' therefore optimal for the maximisation of the Theta values and will be used 
#' by the EM-algorithm. For the multi_HMM_EM(), small changes where made to 
#' calculate the index of the Theta alues right.
#' 
#' @param m number of likelihoods
#' @param N length of the supplied dataset
#' @param u output matrix of the u_function
#' @param x a sample of a Hidden Markov Model 
#' @param theta Theta vector
#' @param L1 likelihood of the first hidden state
#' @param L2 likelihood of the second hidden state
#' @param L3 optional. likelihoods of the third hidden state
#' @param L4 optional. likelihoods of the 4th hidden state
#' @param L5 optional. likelihoods of the 5th hidden state
#' @param start_index index paramter to assign the right amount of thetas 
#' to the likelihoods 
#'
#'@return returns the corresponding likelihood
#'
#'@details For more detailed explanation we recommend the source Hidden Markov 
#'Models for Times Series by Walter Zucchini, Iain MacDonald & Roland Langrock, 
#'especially page 72.



multi_theta_estimation <- function( m, N, u, x, theta, L1, L2, L3 = NULL, 
                                    L4 = NULL,L5 = NULL, start_index ){
  
  #For each state, the individual likelihood is logarithmized and summed with 
  #the coresponding u-vector over all timepoints. The contributions of the 
  #states are then added together to a global likelihood that needs to be 
  #maximized 
  
  #Due to the fact that the maximization is dependent on the number of likelihoods
  #we query for the given number of states 
  
  #The index start_index is used to rightly identify the first and the last 
  #theta-value of the grouping vector Theta (see description in multi_HMM_EM())
  
  l3<-0
  
  
  
  #For the first likelihood
  l3 <- l3 + u[, 1] %*% (log(L1(x, theta[start_index[1]:(start_index[2]-1)])))
  #For the second likelihood
  l3 <- l3 +u[, 2] %*% (log(L2(x, theta[start_index[2]:(start_index[3]-1)])))
  
  #If we have more than two likelihoods we can use this algorithm
  if (m>2){
    l3 <- l3 +u[, 3 ]%*% (log(L3(x, theta[start_index[3]:(start_index[4]-1)])))
  }
  if (m>3) {
    l3 <- l3 +u[, 4] %*% (log(L4(x, theta[start_index[4]:(start_index[5]-1)])))
  }
  if (m>4){
    l3 <- l3 +u[, 5] %*% (log(L5(x, theta[start_index[5]:(start_index[6]-1)])))
  }
  return(l3*-1)
}


