#' Fitting a Hidden Markov Model via Expectation-Maximisation algorithm
#'
#' @description Estimation of the transition probabilites, the initial state probabilites and the hidden state parameters of a Hidden Markov Model
#' by using the Baum-Welch Algorithm.
#'
#' @param x, a sample of a Mixed Model
#' @param m, the number of states
#' @param L1, likelihood of the first hidden state
#' @param L2, likelihood of the second hidden state
#' @param L3-L5, optional. likelihoods of the third, 4th and 5th hidden state
#' @param iterations, optional. Number of iterations for the EM-Algorithm
#' @param delta, optional. Stop criterion for the EM-Algorithm
#'
#' @return The estimated parameters are rounded by 3 decimals and returned in a list.
#' @details This function estimates the hidden states of the Hidden Markov Model with the
#' help of the Baum Welch algorithm. The function iterativly applies both estimation and
#' maximisation steps to arrive at the predicted parameters. When the maximum difference between
#' the functions is abitrarily small (see delta) the iteration stops.
#'





HMM2<-function(x, m, L1, L2, L3=NULL, L4=NULL, L5=NULL, iterations=NULL, delta=NULL){

#setting starting values:

theta_hat<-c()
theta_hat[1:m]<-sample(x, size=m)

Gamma_hat<-matrix(, ncol = m, nrow = m)
Gamma_hat[,]<-1/m

sigma_hat<-matrix(, nrow = m)
sigma_hat[,]<-1/m


#sample length
N<-length(x)

#likelihoods
p1<-L1(x, theta_hat[1])
p2<-L2(x, theta_hat[2])
p<-c(p1,p2)

if(!is.null(L3)){
  p3<-L3(x, theta_hat[3])
  p<-c(p,p3)
}
if(!is.null(L4)){
  p4<-L4(x, theta_hat[4])
  p<-c(p,p4)
}
if(!is.null(L5)){
  p5<-L5(x, theta_hat[5])
  p<-c(p,p5)
}


#We define a treshold of the length of x for every likelihood
set<-seq(1, length(p)-N+1, length.out = m)


#while parameters
if(is.null(iterations)){
  z<-500
} else {
  z<-iterations
}

if(is.null(delta)){
  d<-0.0001
} else {
  d<-delta
}

q<-1
delta<-Inf


while (delta>d && q<z){

  q<-q+1

  p1<-L1(x, theta_hat[1])
  p2<-L2(x, theta_hat[2])
  p<-c(p1,p2)

  if(!is.null(L3)){
    p3<-L3(x, theta_hat[3])
    p<-c(p,p3)
  }
  if(!is.null(L4)){
    p4<-L4(x, theta_hat[4])
    p<-c(p,p4)
  }
  if(!is.null(L5)){
    p5<-L5(x, theta_hat[5])
    p<-c(p,p5)
  }

  theta<-theta_hat
  Gamma<-Gamma_hat
  sigma<-sigma_hat


  out<-alpha_function(m, N, sigma, Gamma, p, set)
  alpha <- out[[1]]
  weight <- out[[2]]

  beta<-beta_function(m, N, Gamma, p, set, c= weight )


  u<-u_function(m, N, alpha, beta)

  v<-v_function(m, N,beta, p,c=weight, Gamma, set, u = u)


  #set in the likelihood (consisting of three parts)
  sigma_hat<-sigma_hat<-u[1,]

  f<-apply(v,2,sum)
  f<-matrix(f, ncol = m, byrow = N)

  Gamma_hat<-matrix(, ncol = m, nrow = m)

  for (j in 1:m){
    Gamma_hat[j,]<-f[j,]/sum(f[j,])
  }

  #We maximize the third part of the likelihood with nlminb
  #Depending on the number of likelihoods

  if(m==2){
    theta_hat<-nlminb(start = theta, theta_estimation, m=m, N=N, u=u, x=x, L1=L1, L2=L2)$par
  } else if (m==3){
    theta_hat<-nlminb(start = theta, theta_estimation, m=m, N=N, u=u, x=x, L1=L1, L2=L2, L3=L3)$par
  } else if (m==4){
    theta_hat<-nlminb(start = theta, theta_estimation, m=m, N=N, u=u, x=x, L1=L1, L2=L2, L3=L3, L4=L4)$par
  } else if (m==5){
    theta_hat<-nlminb(start = theta, theta_estimation, m=m, N=N, u=u, x=x, L1=L1, L2=L2, L3=L3, L4=L4, L5=L5)$par

  }

  #Delta Estimation
  delta1<-max(abs(sigma_hat-sigma))
  delta2<-max(abs(Gamma_hat-Gamma))
  delta3<-max(abs(theta_hat-theta))
  delta<-max(delta1, delta2, delta3)

}

return<-list( "method of estimation:" = "Maximisation via the EM-Algorithm",
              Sigma=round(sigma_hat,3),
              Gamma=round(Gamma_hat,3),
              Theta=round(theta_hat, 3),
              Iterations=q,
              Delta=delta )

return(return)
}

