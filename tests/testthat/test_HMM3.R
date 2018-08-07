context("HMM3 Test control")
source("Transformation.R")

test_that("Probit Control mit p=0.5", {
  #Due to the fact that the probit transformation uses the exponential function
  #we should get the probability 0.5 for Sigma and Gamma alike if we use m=2 
  #and the factor c(0,0,0,1,1). Theta will be a vector of c(1,1)
  
  prob_test <- trans(c(0,0,0,1,1),m=2)
  
  expect_equal(prob_test[,1],c(0.5,0.5))
  expect_equal(prob_test[,2] ,c(0.5,0.5))
  expect_equal(prob_test[,3] ,c(0.5,0.5))
  expect_equal(prob_test[,4],c(1,1))
})




