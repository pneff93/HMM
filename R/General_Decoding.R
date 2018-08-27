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
#' 



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
  
  
  return(list(Local_Decoding = local_path_out,
              Global_Decoding = global_path_out,
              "Differences in decoding" = question))
  }  
  
  
  
