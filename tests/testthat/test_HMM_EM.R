context("single_HMM_EM Test control")

test_that("Output of the Alpha function", {
  #We have to states m = 2 and two x1, x2 thus N = 2
  #Delta and Gamma have all propabilities equal to 0.5 
  #both states are one time very likely (identity matrix) and one time very 
  #unlikely. We expect that the weighted alpha should be equal to the enitial
  #Gamma matrix 
  
  delta_test <- matrix(c(0.5, 0.5), nrow = 2)
  gamma_test <- matrix (c(0.5, 0.5, 0.5, 0.5), 2, 2) 
  p_test <- c(1, 1, 1, 1)
  set_test <- c(1, 3)
  
  alpha_test_1 <- alpha_function(m = 2, N = 2, delta_test, gamma_test, p_test,
                                 set_test)
  
  expect_equal(alpha_test_1[[1]], matrix (c(0.5, 0.5, 0.5, 0.5), 2, 2))
  expect_equal(alpha_test_1[[2]] , c(1, 1))
  
  #Now for the fact that they are fully unlikely, we expect NaN due to the fact,
  #that we devide through zero for the weights the first sum is equal to zero
  #(sums of zero prob.), and the second one should also be NaN.
  p_test_2 <- c(0, 0, 0, 0)
  alpha_test_2 <- alpha_function(m = 2, N = 2, delta_test, gamma_test, p_test_2,
                                 set_test)
  
  expect_equal(alpha_test_2[[1]], matrix (c(NaN, NaN, NaN, NaN), 2, 2))
  expect_equal(alpha_test_2[[2]] , c(0, NaN))
  
})


test_that("Output of Beta function", {
  #Same parameters as the alpha function. The weights are set to 1 
  delta_test <- matrix(c(0.5, 0.5), nrow = 2)
  gamma_test <- matrix (c(0.5, 0.5, 0.5, 0.5), 2, 2) 
  p_test <- c(1, 1, 1, 1)
  set_test <- c(1, 3)
  c_test <- c(1, 1)
  
  beta_test_1 <- beta_function(2, 2, gamma_test, p_test, set_test, c_test)
  
  expect_equal(beta_test_1, matrix(c(1, 1, 1, 1), 2, 2))
})

test_that("Output u-function and v-function", {
  #we set the previous alphas and betas and expect the output to be
  #0.5, 0.5
  alpha_test <- matrix (c(0.5, 0.5, 0.5, 0.5), 2, 2) 
  beta_test <- matrix(c(1, 1, 1, 1), 2, 2)
  gamma_test <- matrix (c(0.5, 0.5, 0.5, 0.5), 2, 2) 
  c_test <- c(1, 1)
  p_test <- c(1, 1, 1, 1)
  set_test <- c(1, 3)
  
  u_test <- u_function(2, 2, alpha_test, beta_test)
  v_test <- v_function(2, 2, beta_test, p_test, c_test, gamma_test, set_test,
                       u = u_test)
  
  
  expect_equal(u_test, matrix (c(0.5, 0.5, 0.5, 0.5), 2, 2) )
  expect_equal(v_test, matrix(c(0.25, 0.25, 0.25, 0.25), 1, 4))
  })


  
  
