context("HMM2 Test control")

test_that("Testing the output of the alpha function", {
  #We have to states m=2 and two x1,x2 thus N=2
  #sigma and gamma have all propabilities equal to 0.5 
  #both states are one time very likely (identity matrix) and one time very unlikely
  #We expect that the weighted alpha should be equal to the enitial gamma matrix 
  sigma1 <-matrix(c(0.5,0.5), nrow = 2)
  gamma1 <- matrix (c(0.5,0.5,0.5,0.5),2,2) 
  p1 <- c(1,1,1,1)
  set1 <- c(1,3)
  
  alpha_test1 <- alpha_function(m=2,N=2,sigma1,gamma1,p1,set1)
  
  expect_equal(alpha_test1[[1]],matrix (c(0.5,0.5,0.5,0.5),2,2))
  expect_equal(alpha_test1[[2]] ,c(1,1))
  
  #Now for the fact that they are fully unlikely, we expect NaN due to the fact, that we devide through zeor 
  #for the weights the first sum is equal to zero (sums of zero prob.), and the second one should also be NaN
  p2 <- c(0,0,0,0)
  alpha_test2 <- alpha_function(m=2,N=2,sigma1,gamma1,p2,set1)
  
  expect_equal(alpha_test2[[1]],matrix (c(NaN,NaN,NaN,NaN),2,2))
  expect_equal(alpha_test2[[2]] ,c(0,NaN))
  
})


test_that("Output of Beta function",{
  #Same parameters as the alpha function. The weights are set to 1 
  sigma1 <-matrix(c(0.5,0.5), nrow = 2)
  gamma1 <- matrix (c(0.5,0.5,0.5,0.5),2,2) 
  p1 <- c(1,1,1,1)
  set1 <- c(1,3)
  
  beta1_test <- beta_function(2,2,gamma1,p1,set1,c=c(1,1))
  
  expect_equal(beta1_test,matrix(c(1,1,1,1),2,2))
})
