context("HMM3 Test control")

test_that("Probit Control mit p = 0.5", {
  #Due to the fact that the probit transformation uses the exponential function
  #we should get the probability 0.5 for Delta and Gamma alike if we use m = 2 
  #and the factor c(0, 0, 0, 1, 1). Theta will be a vector of c(1, 1).
  
  prob_test <- trans(c(0, 0, 0, 1, 1), m = 2)
  
  expect_equal(prob_test[[1]], c(0.5, 0.5))
  expect_equal(prob_test[[2]], matrix(c(0.5, 0.5, 0.5, 0.5), 2, 2))
  expect_equal(prob_test[[3]], c(1, 1))
})

test_that("Control of Likelihood function", {
  #With a fixed set of parameters we expect to gain the same result for the 
  #likelihood every time. factor = c(0, 0, 0, 1, 1),  L1 = L2 = normal distr. 
  #with sd = 1 and m = 2, x1 = c(1, 1, 0, 0) and x2 =c(1, 1, 1, 1)
  expect_equal(round(LH(c(0, 0, 0, 1, 1), m = 2, x = c(1, 1, 0, 0), L1 = dnorm,
                        L2 = dnorm), 4), 4.6758)
  expect_equal(round(LH(c(0,0,0,1,1),m = 2,x = c(1,1,1,1),L1 = dnorm,L2 = dnorm)
                     ,4),3.6758)
})


