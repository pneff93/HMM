
################################################################################
############################## General_Alpha ###################################


#' Normalized Forward Probability function EM
#'
#' @description This function calculates the normalized forward probability of 
#' the HMM states. The probabilities are implemented by calculating first the 
#' forward probabilities of each step and then they are weighted such that the 
#' sum of probabilities are equal to one. This procedure is done to prevent the 
#' threat of underflow. In the general literatur the forward probabilities are
#' also revered to as alpha matrix.
#'
#' @param m number of likelihoods
#' @param N length of the supplied dataset
#' @param delta Delta vector
#' @param Gamma Gamma matrix
#' @param p vector of likelihood probabilities of the dataset
#' @param set index vector to align the p vector
#'
#' @return Returns a list of parameters. The first list are the normalised 
#' alphas as a matrix with the dimension m as columns and the datapoints t as 
#' rows. The second part of the list are the corresponding weights.



alpha_function <- function( m, N, delta, Gamma, p, set ){
  
  alpha<-matrix(, ncol = m, nrow = N)
  weight <- c()
  weight[1] <- sum(t(delta) %*% diag(c(p[set])))
  alpha[1, ]<-(t(delta) %*% diag(c(p[set]))) / weight[1]
  
  for (t in 2:N){
    prod <- alpha[t-1, ] %*% Gamma %*% diag(c(p[set+t-1]))
    weight[t] <- sum(prod)
    alpha [t, ] <- prod / weight[t]
  }
  out <- list(alpha, weight)
  return(out)
}





################################################################################
############################## General_Beta ####################################


#' Normalized Backward Probability function EM
#'
#' @description This function calculates the normalized backward probability 
#' of the HMM states. For this the function needs the weights of the forward
#' probability, which have to be inputed as parameter weight. 
#' In the general literatur also revered to as beta matrix.
#'
#' @param m number of likelihoods
#' @param N length of the supplied dataset
#' @param Gamma Gamma matrix
#' @param p vector of likelihood probabilities of the dataset
#' @param set index vector to align the p vector
#' @param weight weights provided by the alpha calculation
#'
#' @return Returns the beta as a matrix with the dimensions m as rows
#' and the datapoints t as columns.



beta_function <- function( m, N, Gamma, p, set, weight ){
  
  beta<-matrix(, ncol = N, nrow = m)
  beta[, N]<-rep(1, times = m) / weight[N]
  
  for (t in (N-1):1){
    beta[, t] <- Gamma %*% diag(c(p[set+t])) %*% beta[, t+1] / weight[t]
  }
  
  return(beta)
}





################################################################################
############################## General_Vector_u ################################


#' Normalized u-function EM
#'
#' @description This function calculates the u-function which is part of the 
#' estimation step. Each element in u aligns with the conditional expectation to
#' reach state j given your datapoint x(t). The corresponding weights are 
#' already included in the alpha and beta values.
#'
#'
#' @param m number of likelihoods
#' @param N length of the supplied dataset
#' @param alpha normalized alpha matrix
#' @param beta normalized beta matrix
#'
#' @return returns the corresponding normalized u's



u_function <- function( m, N, alpha, beta ){
  
  
  u <- matrix(, ncol = m, nrow = N)
  
  for (t in 1:N){
    for (j in 1:m){
      u[t, j] <- (alpha[t, j] %*% t(beta[j, t])) / (alpha[t, ] %*% beta[, t])
    }
  }
  return(u)
}





################################################################################
############################## General_Vector_v ################################


#' Normalized v-function EM
#'
#' @description This function calculates the v-function which is part of the 
#' estimation step. Each element in v aligns with the conditional expectation 
#' to reach state k from the previous state j given your datapoint x(t).
#' Due to the normalizing of the factors we also need to include the
#' corresponding weights into our function.
#'
#'
#' @param m number of likelihoods
#' @param N length of the supplied dataset
#' @param beta beta matrix
#' @param p vector of likelihood probabilities of dataset
#' @param weight weights provided by the alpha calculation
#' @param Gamma Gamma matrix
#' @param set index vector to align the p vector
#' @param u output matrix of the u_function
#'
#' @return returns the coresponding values, but due to further calculation steps
#' the m x m matrix for each timepoint of v-values is returned as one row in 
#' the output matrix.



v_function <- function( m, N, beta, p, weight, Gamma, set, u ){
  
  
  v <- matrix(, nrow = 0, ncol = m*m)
  out <- matrix(, nrow = m, ncol = m)
  
  for (t in 2:N){
    pp <- diag(c(p[set+t-1]))
    
    for (i in 1:m){
      for(j in 1:m){
        out[i, j] <- ((u[t-1, i] * Gamma[i, j] * pp[j, j] * beta[i, t])
                      / (beta[i, t-1] * weight[t-1]))
        
      }
    }
    v <- c(v, c(t(out)))
    
  }
  v <- matrix(v, nrow = N-1, byrow = N)
  
  return(v)
}





################################################################################
############################## General_Transformation ##########################


#' Transformation function - DM
#' 
#' @description The transformation function transforms a non-restricted parameter
#' vector into a restricted Delta, Gamma and Theta output.  
#' 
#' @param factor vector of unrestricted parameters
#' @param m number of Likelihoods 
#' 
#' @return List of Delta, Gamma, Theta
#' 
#' @details In the direct maximisation the nlminb()- minimisation function can not 
#' directly implement the constraints of the parameter values. 
#' This function ensures that the estimated parameters of the direct 
#' optimisation still fullfill their requirements. The requirements are that all
#' probabilities are between zero and one and that the rows of Gamma (as well as
#' the vector Delta) sum up to one. For this transformation the logit model is 
#' used. 
#' 
#' From the input vector factor, the elements are extracted as follows and used 
#' for transformation: 
#' 
#' Delta vector (1 x m): The logit transformation requires (m-1)  elements, thus 
#' the first (m-1) elements of the factor vector are used. 
#'  
#' Gamma matrix (m x m): The logit transformation requires per row (m-1) elements
#' thus in total m(m-1) elements required. Accordingly the elements from index 
#' m til (m+1)(m-1) are extracted from the factor vector. 
#' 
#' Theta vector: Due to the fact that the trans() function can be used both from 
#' the single and multifactor Direct Maximisation the Theta vector has no pre-
#' defined lenght, but needs to have at least the length of m (thus each 
#' likelihood has at least one parameter). 
#' 
#' The resulting Delta, Gamma and Theta values are returned as a single list. 



trans <- function( factor, m ){
  ##########################
  #Transformation Delta-vector
  
  #Building the constrains for Delta via logit transformation:
  # delta[i] <- exp(factor[i])/sum(1+exp(factor)) for i = 1 - (m-1)
  # delta[m] <- 1/sum(1+exp(factor))
  
  delta <- c()
  div <- 1 + sum(exp(factor[1:(m-1)]))
  for (i in 1:(m-1)){
    delta[i] <- exp(factor[i]) / div
  }
  delta[m] <- 1 / div
  
  
  ##########################
  #Transformation Gamma-matrix
  
  #Building restrictions on Gamma we extract therefore the next m(m-1) elements
  #of factor; so the elements: (m-1)+1 = m till (m-1)+m(m-1) = (m+1)(m-1)
  
  #For the Gamma matrix we determine the off-diagonal elements: 
  
  #for i =! j: 
  #exp(factor[i])/sum(1+exp(factor))
  #for i == j 
  #1/sum(1+exp(factor))
  
  #As example for m = 3 
  # ( --- tau12 tau13)
  # (tau21 ---  tau23)
  # (tau31 tau32 ---)
  
  
  #Recall that each rowsum has to be equal to one 
  #The variable count ensures that only the off-diagonal elements are calculated
  #with the factor vector 
  
  x <- matrix(factor[m:((m+1)*(m-1))], ncol = (m-1), nrow = m, byrow = TRUE)
  
  
  gamma <- matrix(, nrow = m, ncol = m)
  
  for (i in 1:m){
    div <- 1+sum(exp(x[i, ]))
    count <- 0
    
    for (j in 1:m){
      if(i==j){ 
        gamma[i, j] <- 1/div
      } else {
        count <- count + 1
        gamma[i, j] <- exp(x[i, count])/div
      }
    }
  }
  
  ##########################
  #Transformation Theta-values
  
  
  #The Theta values are not manipulated and are just stored in an individual 
  #vector 
  
  theta <- factor[((m-1)*(m+1)+1):length(factor)]
  
  
  return(list(delta, gamma, theta))
  }





################################################################################
############################## General_Decoding ################################


#' Local and global decoding 
#' 
#' @description Calculation of local and global decoding given the dataset and
#' the coresponding parameters of the HMM
#' 
#' 
#' @param x a sample of a Hidden Markov Model 
#' @param m the number of states
#' @param L1 likelihood of the first hidden state
#' @param L2 likelihood of the second hidden state
#' @param L3 optional. likelihood of the third hidden state
#' @param L4 optional. likelihood of the 4th hidden state
#' @param L5 optional. likelihood of the 5th hidden state
#' @param gamma estimated Gamma matrix from previous estimation
#' @param delta estimated Delta vector from previous estimation
#' @param theta estimated Theta values from previous estimation
#' @param multi parameter. TRUE if multifactor HMM() is used 
#' 
#' @return Returns the path of the local and global decoding as well as a
#' statement, whether the two paths differ.



decode <- function( x, m, L1, L2, L3 = NULL, L4 = NULL, L5 = NULL, gamma,
                    delta, theta, multi = FALSE ){
  
  ##########################
  #conducting the probability vector 
  
  if (multi==FALSE){
    theta_hat <- theta 
    p1 <- L1(x, theta_hat[1])
    p2 <- L2(x, theta_hat[2])
    p <- c(p1, p2)
    
    if(!is.null(L3)){
      p3 <- L3(x, theta_hat[3])
      p <- c(p, p3)
    }
    if(!is.null(L4)){
      p4 <- L4(x, theta_hat[4])
      p <- c(p, p4)
    }
    if(!is.null(L5)){
      p5 <- L5(x, theta_hat[5])
      p <- c(p, p5)
    }
  } else {
    
    theta_hat <- c()
    set2 <- c()
    #rewriting the Thetas (see the multi functions for further explanation)
    for (i in 1:length(theta)){
      set2[i] <- length(theta_hat)+1
      theta_hat <- c(theta_hat, theta[[i]])
    }
    
    
    set2[length(set2)+1] <- length(theta_hat)+1
    
    #p-Vector :
    
    p1 <- L1(x, theta_hat[set2[1]:(set2[2]-1)])
    p2 <- L2(x, theta_hat[set2[2]:(set2[3]-1)])
    p <- c(p1, p2)
    
    if(!is.null(L3)){
      p3 <- L3(x, theta_hat[set2[3]:(set2[4]-1)])
      p <- c(p, p3)
    }
    if(!is.null(L4)){
      p4 <- L4(x, theta_hat[set2[4]:(set2[5]-1)])
      p <- c(p, p4)
    }
    if(!is.null(L5)){
      p5 <- L5(x, theta_hat[set2[5]:(set2[6]-1)])
      p <- c(p, p5)
    }    
  }
  
  
  
  ##########################
  #computing forward and backward probabilities
  #computing the alpha and beta matrices and the Likelihood L_T that 
  #is the sum of the weigths at timepoint T
  
  N <- length(x)
  set <- seq(1, m*N-N+1, length.out = m)
  
  out <- alpha_function( m, N, delta, gamma, p, set )
  alpha <- out[[1]]
  weight <- out[[2]]
  
  beta <- beta_function( m, N, gamma, p, set, weight = weight )
  
  
  #L_T is the sum of the last alpha, which is saved in weight[T]
  L_T <- weight[N]
  
  
  
  ##########################
  #Local decoding: 
  #With Pr(X_t = x_t, C_t = i) = alpha_t(i)* beta_t(i)
  #the result is Pr (C_t = i | X_t  = x_t) =  alpha_t(i)* beta_t(i)/L_T
  
  
  #computation of the mutliplication matrix by multipling each case(i)
  #instead of using a for-loop we use as a trick the simple scalar 
  #multilpication such that the elements are multiplied with each other.
  #This enhances the performance.
  
  probmatrix <- matrix(, nrow = N, ncol = m)
  for (i in 1:m){
    probmatrix[, i] <- alpha[, i] * beta[i, ] /L_T
  }
  
  #for each row we now want the most likely state as return
  #(which is the column name)
  
  local_path_out <- apply(probmatrix, 1, which.max)
  
  
  
  ##########################
  #Global decoding
  
  #Searching for the most likely branch via
  #using the viterbi algorithm:
  
  etha <- matrix(, nrow = N, ncol = m)
  
  #etha_1_i = delta_i * p_i(x1)
  #thus again we use the simple scalar multiplication
  
  #Note, that again we use the normalization to tackle the problem of underflow
  #with the weight w 
  
  #t = 1 
  w <- delta * p[set]
  etha[1, ] <- w / sum(w)
  
  
  
  
  #t = 2, ..., N
  #for higher t we have to calculate the etha by taking the max previous one : 
  #etha_t_j = (max(etha_t-1_i *gamma_ij)*p_j_(x_t))
  #again with scalar multiplication we get by (etha[i, ]*gamma) the corresponding
  #matrix and take the max value of each column, thus selecting the 
  #max (i) for each state j then we multiply each max with the corresponding
  #probability
  
  for ( i in 2:N){
    w <- apply(etha[(i-1), ]*gamma, 2, max) * p[set+i-1]
    etha[i, ] <- w / sum(w)
  }
  
  #after computing the etha we now have to reconstruct the optimal path 
  #recursevly with:
  #i_T = argmax etha_Ti
  #i_t = argmax (etha_ti * gamma_i, i+1)
  
  #(again with multiplying as a scalar)
  
  global_path_out <- vector(, length = N)
  global_path_out[N] <- which.max(etha[N, ])
  
  for ( i in (N-1):1){
    global_path_out[i] <- which.max(gamma[, global_path_out[i+1]]*etha[i, ])
  }
  
  
  
  question <- (sum(abs(global_path_out-local_path_out)) != 0)
  
  output <- list(Local_Decoding = local_path_out,
                 Global_Decoding = global_path_out,
                 "Differences in decoding" = question)
  
  return(output)
}  


