##RESULTS:
#Both methods estimate the parameters with high speed and accurately



#transition matrix
gamma<-matrix(c(0.8,0.2,0.3,0.7), byrow=T, nrow=2)

#initial state probabilities
delta<-c(0.3,0.7)

#sample size
n<- 1000

x<-c()
set.seed(10)
s1<-rpois(10000,10)
s2<-rpois(10000, 2)

#initial state
random_number<-runif(1, 0, 1)

if (random_number < delta[1]){
  p<-1
  x[1]<-sample(s1, 1, replace = F)
} else {
  p<-2
  x[1]<-sample(s2, 1, replace = F)
}

#sample creation
for (i in 2:n){
  random_number<-runif(1, 0, 1)
  if (random_number < gamma[p,1]){
    p<-1
    x[i]<-sample(s1, 1, replace = F)
  } else {
    p<-2
    x[i]<-sample(s2, 1, replace = F)
  }
}


#likelihoods
L1<-function(x, lambda){
  p1<-lambda^x/factorial(x)*exp(-lambda)
  return(p1)
}


L2<-function(x, lambda){
  p2<-lambda^x/factorial(x)*exp(-lambda)
  return(p2)
}

m<-2

HMM(x=x,m=m,method="EM",L1=L1,L2=L2)
HMM(x=x,m=m,method="DM",L1=L1,L2=L2)

