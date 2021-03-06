context("multi_HMM_DM Test control")

#This test section is for controlling the multi_HMM_DM function. 
#We do not control the trans-function, because this is already taken 
#care of in the HMM3-test.


#First test is about the multi_LH function 

test_that("multi_LH function control", {
  
    #We use the same likelihoods as the mulitHMM2 test function.
    #with a fixed set of parameters we expect to gain the same result for the 
    #likelihood every time. factor = c(0, 0, 0, 1, 1), L1 = L2 = normal distr. 
    #with mu = 1, sd = 1 and m = 2, x1 = c(1, 1, 0, 0) and x2 = c(1, 1, 1, 1)
    #the result should be the similar output like the LH-test from HMM3-test.
    
    L1 <- function(x, theta){
      mu <- theta[1]
      sd <- theta[2]
      p1 <- 1/sqrt(2*pi*(sd^2)) * exp(-((x-mu)^2) / (2*sd^2))
      return(p1)
    }
    
    L2 <- function(x,theta){
      mu <- theta[1]
      sd <- theta[2]
      p2 <- 1/sqrt(2*pi*(sd^2)) * exp(-((x-mu)^2) / (2*sd^2))
      return(p2)
    }
    
    factor_test_1 <- c(0, 0, 0, 1, 1, 1, 1)
    factor_test_2 <- c(0, 0, 0, 1, 1, 1, 1)
    start_index_1 <- c(1, 3, 5) 
    
    expect_equal(round(multi_LH(factor = factor_test_1, m = 2, x = c(1, 1, 0, 0),
                                L1, L2, start_index = start_index_1), 4), 4.6758)
    expect_equal(round(multi_LH(factor = factor_test_2, m = 2, x = c(1, 1, 1, 1),
                                L1, L2, start_index = start_index_1), 4), 3.6758)
    
    
    
})
  

  #Testing for the wrong input in the multi_HMM_DM function, thus if not list is 
  #provided same test format like the multi_HMM_EM test file 

test_that("Differents multiple theta inputs", {
    z <- c(0, 0, 1, 1)
    
    L1 <- function(x, theta){
      mu <- theta[1]
      sd <- theta[2]
      p1 <- 1/sqrt(2*pi*(sd^2)) * exp(-((x-mu)^2) / (2*sd^2))
      return(p1)
    }
    
    L2 <- function(x,theta){
      mu <- theta[1]
      sd <- theta[2]
      p2 <- 1/sqrt(2*pi*(sd^2)) * exp(-((x-mu)^2) / (2*sd^2))
      return(p2)
    }
    
    #theta_test_1 is the right starting input for our estimation
    theta_test_1 <- list(c(0, 1), c(0, 1))
    
    #theta_test_2 is not a list, thus should return an warning
    theta_test_2 <- c(0, 1) 
    
    #theta_test_3 is a list but has not the right amount of parameters, 
    #thus should retur an warning
    theta_test_3 <- list(1, c(0, 1))
    
    output_1 <- multi_HMM_DM(x = z, theta_test_1, m = 2, L1, L2)
    

    
    #theta_test_1
    #The estimated Thetas should be the same and have a mu = 0.5 and sigma = 0.5 
    expect_equal(output_1$Theta[[1]], c(0.5, 0.5))
    expect_equal(output_1$Theta[[2]], c(0.5, 0.5))
    
    
    #compared to the EM-algorithm our code does not produce an error, because 
    #the nlminb() function does not find the right parameters and just returns
    #the input parameters. BUT it does not produce an error, only a warning.
    
    
    #theta_test_2
    expect_warning(multi_HMM_DM(x = z, theta_test_2, m = 2, L1, L2)) 
    #theta_test_3 
    expect_warning(multi_HMM_DM(x = z, theta_test_3, m = 2, L1, L2))
})