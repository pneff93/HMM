##RESULTS:
#Both methods estimate the parameters with high speed and accurately



#transition matrix
gamma <- matrix(c(0.8, 0.2, 0.3, 0.7), byrow = TRUE, nrow = 2)

#initial state probabilities
delta <- c(0.3, 0.7)

#sample size
n <- 1000

x <- c()
set.seed(10)
s1 <- rnorm(10000, 7, 1)
s2 <- rnorm(10000, 2, 10)

#initial state
random_number <- runif(1, 0, 1)

if (random_number < delta[1]){
  p <- 1
  x[1] <- sample(s1, 1, replace = FALSE)
} else {
  p <- 2
  x[1] <- sample(s2, 1, replace = FALSE)
}

#sample creation
for (i in 2:n){
  random_number <- runif(1, 0, 1)
  if (random_number < gamma[p,1]){
    p <- 1
    x[i] <- sample(s1, 1, replace = FALSE)
  } else {
    p <- 2
    x[i] <- sample(s2, 1, replace = FALSE)
  }
}


#likelihoods
L1 <- function(x, mu){
  p1 <- 1/sqrt(2*pi) * exp(-0.5*(x-mu)^2)
  return(p1)
}

L2 <- function(x, theta){
  mu <- theta[1]
  sd <- theta[2]
  p2 <- 1/sqrt(2*pi*(sd^2)) * exp(-((x-mu)^2) / (2*sd^2))
  return(p2)
}

theta1 <- list(10, c(0, 5))
theta1
m <- 2


multi_HMM(x = x, theta = theta1, m = m, method = "EM",
               L1 = L1, L2 = L2, decoding = TRUE)
multi_HMM(x = x, theta = theta1, m = m, method = "DM",
               L1 = L1, L2 = L2, decoding = TRUE)

