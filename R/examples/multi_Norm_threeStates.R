################
#First: Generating the sample of the HMM with the true following true values:

#transition matrix
gamma <- matrix(c(0.7, 0.25, 0.05,
                0.1, 0.4, 0.5,
                0.3, 0.5, 0.2), byrow = TRUE, nrow = 3)

#initial state probabilities
delta <- c(0.5, 0.3, 0.2)

#sample size
n <- 500

#number of states 
m <- 3

#sampling from normal distribution with different mu's and sigma's: 
x <- c()
set.seed(100)
s1 <- rnorm(10000, 7, 1)
s2 <- rnorm(10000, 2, 4)
s3 <- rnorm(10000, 15, 3)


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

#To show complexity the first likelihood only requires one parameter,
#while the other two need two. Theta is always a vector. 

L1 <- function(x, theta){
  mu <- theta 
  p1 <- 1/sqrt(2*pi) * exp(-0.5*(x-mu)^2)
  return(p1)
}

L2 <- function(x, theta){
  mu <- theta[1]
  sd <- theta[2]
  p2 <- 1/sqrt(2*pi*(sd^2)) * exp(-((x-mu)^2) / (2*sd^2))
  return(p2)
}

L3 <- function(x, theta){
  mu <- theta[1]
  sd <- theta[2]
  p3 <- 1/sqrt(2*pi*(sd^2)) * exp(-((x-mu)^2) / (2*sd^2))
  return(p3)
}
################
#Third: Guessing the initial Theta and execution of multifactor function
#intial estimates of Theta
theta1 <- list(8, c(1, 1), c(20, 1))

#execution of both multifactorHMM functions, with decoding=TRUE 
multifactorHMM(x = x, theta = theta1, m = m, method = "EM",
               L1 = L1, L2 = L2, L3 = L3, decoding = TRUE)
multifactorHMM(x = x, theta = theta1, m = m, method = "DM",
               L1 = L1, L2 = L2, L3 = L3, decoding = TRUE)

