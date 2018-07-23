#' Fitting a Hidden Markov Model
#' @description Estimation of the transition probabilites, the initial state probabilites and the hidden state parameters of a Hidden Markov Model
#' by using an EM-Algorithm
#' @param x, a sample of a Mixed Model
#' @param m, the number of states
#' @param L1, likelihood of the first hidden state
#' @param L2, likelihood of the second hidden state
#' @param L3-L5, optional. likelihoods of the third, 4th and 5th hidden state
#' @param iterations, optional. Number of iterations for the EM-Algorithm
#' @param delta, optional. Stop criterion for the EM-Algorithm
#' @return The estimated parameters.
#' @details Using an EM-Algorithm we have bla formula, backwards, forwards
#' @export
#'
HMM<-function(x, m, L1, L2, L3=NULL, L4=NULL, L5=NULL, iterations=NULL, delta=NULL){


#Idee statt direkt einer gro?en Diagonalmatrix mit vielen Zeros
#speichern wir es als vector und ziehen uns dann immer die passenden Elemente
#und basteln nur mit diesen eine (kleine) Diagonalmatrix
#Dadurch Anzahl der states egal !!!

#Initial values

theta_hat<-c()
theta_hat[1:m]<-sample(x, size=m)

Gamma_hat<-matrix(, ncol = m, nrow = m)
Gamma_hat[,]<-1/m

sigma_hat<-matrix(, nrow = m)
sigma_hat[,]<-1/m


#L?nge des samples
N<-length(x)

#Likelihoods
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


#wo beginnt neues p[i]
set<-seq(1, length(p)-N+1, length.out = m)

#f?r s= , ziehe ich mir dann die Werte und bastel eine Diagonalmatrix
#diag(c(p[set+s-1]))


#while parameters
if(is.null(iterations)){
  z<-500
} else {
  z<-iterations
}

if(is.null(delta)){
  d<-0.01
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


  alpha<-alpha_function( m, N, sigma, Gamma, p, set)
  beta<-beta_function( m, N, Gamma, p, set )

  L<-Likelihood(alpha, beta)

  u<-u_function(m, N, alpha, beta, L)
  v<-v_function(m, N, alpha, beta, p, L, Gamma, set)

  #Einsetzen in die dreier Likelihood solutions
  sigma_hat<-sigma_hat<-u[1,]

  f<-apply(v,2,sum)
  f<-matrix(f, ncol = m, byrow = N)

  Gamma_hat<-matrix(, ncol = m, nrow = m)

  for (j in 1:m){
    Gamma_hat[j,]<-f[j,]/sum(f[j,])
  }


  #dritte Likelihood

  if(m==2){
    theta_hat<-nlminb(start = theta, theta_estimation, m=m, N=N, u=u, x=x, L1=L1, L2=L2)$par
  } else if (m==3){
    theta_hat<-nlminb(start = theta, theta_estimation, m=m, N=N, u=u, x=x, L1=L1, L2=L2, L3=L3)$par
  } else if (m==4){
    theta_hat<-nlminb(start = theta, theta_estimation, m=m, N=N, u=u, x=x, L1=L1, L2=L2, L3=L3, L4=L4)$par
  } else if (m==5){
    theta_hat<-nlminb(start = theta, theta_estimation, m=m, N=N, u=u, x=x, L1=L1, L2=L2, L3=L3, L4=L4, L5=L5)$par

  }

  #Delta Berechnung
  delta1<-max(abs(sigma_hat-sigma))
  delta2<-max(abs(Gamma_hat-Gamma))
  delta3<-max(abs(theta_hat-theta))
  delta<-max(delta1, delta2, delta3)

}
AIC<--2*log(L)+2*(m+m*m+m)

return<-list(sigma=round(sigma_hat,2), Gamma=round(Gamma_hat,2), theta=round(theta_hat, 2), iterations=q, delta=delta, AIC=AIC)
return(return)
}
