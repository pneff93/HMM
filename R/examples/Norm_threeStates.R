################
#First: Generating the sample of the HMM with the true following true values:

#transition matrix
gamma <- matrix(c(0.9, 0.05, 0.05,
                0.1, 0.4, 0.5,
                0.3, 0.5, 0.2), byrow = TRUE, nrow = 3)

#initial state probabilities
delta <- c(0.5, 0.3, 0.2)

#sample size
n <- 500


#number of states 
m <- 3

#sampling from normal distribution with different mu's but same sigma's:
x <- c()
set.seed(100)
s1 <- rnorm(10000, 7, 1)
s2 <- rnorm(10000, 2, 1)
s3 <- rnorm(10000, 12, 1)


#initial state
random_number <- runif(1, 0, 1)

if (random_number < delta[1]){
  x[1] <- sample(s1, 1, replace = FALSE)
  p <- 1
} else if (random_number < sum(delta[1:2]) && random_number > delta[1]) {
  x[1] <- sample(s2, 1, replace = FALSE)
  p <- 2
} else {
  x[1] <- sample(s3, 1, replace = FALSE)
  p <- 3
}

#sample creation
for (i in 2:n){
  random_number <- runif(1, 0, 1)
  if (random_number < gamma[p,1]){
    p <- 1
    x[i] <- sample(s1, 1, replace = FALSE)
  } else if(random_number < sum(gamma[p,1:2]) && random_number > gamma[p,1]) {
    p <- 2
    x[i] <- sample(s2, 1, replace = FALSE)
  } else {
    p <- 3
    x[i] <- sample(s3, 1, replace = FALSE)
  }
}

#Display of the sample
hist(x)

################
#Second: Defining the likelihoods.


#likelihoods
L1 <- function(x, mu){
  p1 <- 1/sqrt(2*pi) * exp(-0.5*(x-mu)^2)
  return(p1)
}

L2 <- function(x, mu){
  p2 <- 1/sqrt(2*pi) * exp(-0.5*(x-mu)^2)
  return(p2)
}

L3 <- function(x, mu){
  p3 <- 1/sqrt(2*pi) * exp(-0.5*(x-mu)^2)
  return(p3)
}

################
#Third: Executing the two HMM functions

HMM(x = x, m = m, method = "EM", L1 = L1, L2 = L2, L3 = L3, decoding = TRUE)

HMM(x = x, m = m, method = "DM", L1 = L1, L2 = L2, L3 = L3, decoding = TRUE)

