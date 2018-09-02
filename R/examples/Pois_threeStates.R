##RESULTS:
#The "EM" method needs a very long estimation time
#but the results are better compared to the "DM" method


#transition matrix
gamma <- matrix(c(0.4, 0.2, 0.4,
                0.1, 0.4, 0.5,
                0.3, 0.4, 0.3), byrow = TRUE, nrow = 3)

#initial state probabilities
delta <- c(0.5, 0.3, 0.2)

#sample size
n <- 500

x <- c()
set.seed(100)
s1 <- rpois(10000, 10)
s2 <- rpois(10000, 2)
s3 <- rpois(10000, 4)


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


#likelihoods
L1 <- function(x, lambda){
  p1 <- lambda^x / factorial(x) * exp(-lambda)
  return(p1)
}

L2 <- function(x, lambda){
  p2 <- lambda^x / factorial(x) * exp(-lambda)
  return(p2)
}

L3 <- function(x, lambda){
  p3 <- lambda^x / factorial(x) * exp(-lambda)
  return(p3)
}

hist(x)
#number of states 
m <- 3

HMM(x = x, m = m, method = "EM", L1 = L1, L2 = L2, L3 = L3)

HMM(x = x, m = m, method = "DM", L1 = L1, L2 = L2, L3 = L3)

