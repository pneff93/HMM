context("Test of decoding")

test_that("Test if the decoding estimates the right proberties", {
  
  z <- c(0, 0, 10, 10, 0, 10)
  
  L1 <- function(x, mu){
    p1 <- 1/sqrt(2*pi) * exp(-0.5*(x-mu)^2)
    return(p1)
  }
  
  #Test input
  gamma_test <- matrix (c(0.5, 0.5, 0.5, 0.5), 2, 2)
  delta_test <- c(1, 0)
  theta_test <- c(0, 10)
  
  #We expect both paths to return an path of 1, 1, 2, 2, 1, 2, because with Delta
  #c(1, 0) we set the first Theta as first path and the
  #z falls exactly onto the estimated Theta (thus highest probabilty)
  #The decoding styles should not return different results 
  
  test_out <- decode(x = z, m = 2, L1 = L1, L2 = L1, gamma =  gamma_test,
                     delta = delta_test, theta = theta_test, multi = FALSE)
    
  expect_equal(test_out$Local_Decoding, c(1, 1, 2, 2, 1, 2))
  expect_equal(test_out$Global_Decoding, c(1, 1, 2, 2, 1, 2))
  expect_equal(test_out$Local_Decoding, test_out$Global_Decoding)
  
  })