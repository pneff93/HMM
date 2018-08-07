##RESULTS:
#Both methods need a long time to estimate the parameters
#However the results are not satisfying


#transition matrix
gamma<-matrix(c(0.2,0.1,0.4,0.15,0.15,
                0.1,0.1,0.01,0.5,0.29,
                0.3,0.3,0.05,0.3,0.05,
                0.2,0.05,0.35,0.05,0.35,
                0.6,0.01,0.01,0.08,0.3), byrow=T, nrow=5)

#initial state probabilities
sigma<-c(0.5, 0.1, 0.2, 0.1, 0.1)

#sample size
n<- 1000

x<-c()
set.seed(100)
s1<-rpois(10000,10)
s2<-rpois(10000, 2)
s3<-rpois(10000, 4)
s4<-rpois(10000, 8)
s5<-rpois(10000, 1)

#initial state
random_number<-runif(1, 0, 1)

if (random_number < sigma[1]){
  x[1]<-sample(s1, 1, replace = F)
  p<-1
} else if (random_number < sum(sigma[1:2]) && random_number > sigma[1]) {
  x[1]<-sample(s2, 1, replace = F)
  p<-2
} else if (random_number < sum(sigma[1:3]) && random_number > sum(sigma[1:2])){
  x[1]<-sample(s3, 1, replace = F)
  p<-3
} else if (random_number < sum(sigma[1:4]) && random_number > sum(sigma[1:3])){
  x[1]<-sample(s4, 1, replace = F)
  p<-4
} else {
  x[1]<-sample(s5, 1, replace = F)
  p<-5
}

#sample creation
for (i in 2:n){
  random_number<-runif(1, 0, 1)
  if (random_number < gamma[p,1]){
    p<-1
    x[i]<-sample(s1, 1, replace = F)
  } else if(random_number < sum(gamma[p,1:2]) && random_number > gamma[p,1]) {
    p<-2
    x[i]<-sample(s2, 1, replace = F)
  } else if (random_number < sum(gamma[p,1:3]) && random_number > sum(gamma[p,1:2])){
    p<-3
    x[i]<-sample(s3, 1, replace = F)
  } else if (random_number < sum(gamma[p,1:4]) && random_number > sum(gamma[p,1:3])){
    p<-4
    x[i]<-sample(s4, 1, replace = F)
  } else {
    p<-5
    x[i]<-sample(s5, 1, replace = F)
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

L3<-function(x, lambda){
  p3<-lambda^x/factorial(x)*exp(-lambda)
  return(p3)
}

L4<-function(x, lambda){
  p4<-lambda^x/factorial(x)*exp(-lambda)
  return(p4)
}

L5<-function(x, lambda){
  p5<-lambda^x/factorial(x)*exp(-lambda)
  return(p5)
}

m<-5

HMM(x=x,m=m,method="EM",L1=L1,L2=L2,L3=L3,L4=L4,L5=L5)
#HMM(x=x,m=m,method="DM",L1=L1,L2=L2,L3=L3,L4=L4,L5=L5)

