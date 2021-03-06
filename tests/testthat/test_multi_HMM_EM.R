context("multi_HMM_EM Test control")

#Due to the fact that the alpha and beta functions are already tested in the 
#test_HMM2 file this test section focus on the implementation of the multiple
#thetas. For alpha & beta unit test, please refere to test_HMM2.R

#We test the overall multi_HMM_EM function on functinality given different input
#parameters. We have a fixed set of x and two likelihoods (normally-distributed)

test_that("Differents multiple theta inputs", {
  z <- c(0, 0, 1, 1)
  
  L1 <- function(x, theta){
    mu <- theta[1]
    sd <- theta[2]
    p1 <- 1/sqrt(2*pi*(sd^2)) * exp(-((x-mu)^2) / (2*sd^2))
    return(p1)
  }
  
  L2 <- function(x, theta){
    mu <- theta[1]
    sd <- theta[2]
    p2 <- 1/sqrt(2*pi*(sd^2)) * exp(-((x-mu)^2) / (2*sd^2))
    return(p2)
  }
  
  #theta_test_1 is the right starting input for our estimation
  theta_test_1 <- list(c(0, 1), c(0, 1))
  #theta_test_2 is not a list, thus should return an error 
  theta_test_2 <- c(0, 1) 
  #theta_test_3 is a list but has not the right amount of parameters, thus should 
  #retur an error
  theta_test_3 <- list(1, c(0, 1))
  
  output_1 <- multi_HMM_EM(x = z, theta = theta_test_1, m = 2, L1, L2)
  
  #theta_test_1
  #The estimated thetas should be the same and have a mu = 0.5 and sigma = 0.5 
  expect_equal(output_1$Theta[[1]], c(0.5, 0.5))
  expect_equal(output_1$Theta[[2]], c(0.5, 0.5))
  
  
  #theta_test_2
  #note that we suppress the warnings, because the statement causes both error 
  #and warnings 
  expect_error(suppressWarnings(multi_HMM_EM(x = z, theta = theta_test_2, m = 2,
                                             L1, L2))) 
  #theta_test_3
  expect_error(suppressWarnings(multi_HMM_EM(x = z, theta = theta_test_3, m = 2,
                                             L1, L2)))
})
