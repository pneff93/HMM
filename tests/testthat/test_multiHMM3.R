context("multiHMM3 Test control")

#This test section is for controlling the multiHMM3 function. 
#We do not controll the trans-function, because this is already taken 
#care of in the HMM3-test


#First test is about the multiLH function 

test_that("multiLH function control",{
  
    #We use the same likelihoods as the mulitHMM2 test function.
    #with a fixed set of parameters we expect to gain the same result for the likelihood
    #every time. factor= c(0,0,0,1,1), L1=L2=normal distr. with mu=1,sd=1 and m=2 
    # x1 =c(1,1,0,0) and x2 =c(1,1,1,1)
    #the result should be the similar output like the LH-test from HMM3-test
    
    L1<-function(x, theta){
      mu <- theta[1]
      sd <- theta[2]
      p1<-1/sqrt(2*pi*(sd^2))*exp(-((x-mu)^2)/(2*sd^2))
      return(p1)
    }
    
    L2<-function(x,theta){
      mu <- theta[1]
      sd <- theta[2]
      p2<-1/sqrt(2*pi*(sd^2))*exp(-((x-mu)^2)/(2*sd^2))
      return(p2)
    }
    
    factor1 <- c(0,0,0,1,1,1,1)
    factor2 <- c(0,0,0,1,1,1,1)
    set2_1 <-c(1,3,5) 
    
    expect_equal(round(multiLH(factor = factor1,m=2,x=c(1,1,0,0),L1,L2,set2 = set2_1),4), 4.6758)
    expect_equal(round(multiLH(factor = factor2,m=2,x=c(1,1,1,1),L1,L2,set2=set2_1),4), 3.6758)
    
    
    
})
  

  #Testing for the wrong input in the multiHMM3 function, thus if not list is provided
  #same test format like the multiHMM2 test file 

test_that("Differents multiple theta inputs",{
    z <- c(0,0,1,1)
    
    L1<-function(x, theta){
      mu <- theta[1]
      sd <- theta[2]
      p1<-1/sqrt(2*pi*(sd^2))*exp(-((x-mu)^2)/(2*sd^2))
      return(p1)
    }
    
    L2<-function(x,theta){
      mu <- theta[1]
      sd <- theta[2]
      p2<-1/sqrt(2*pi*(sd^2))*exp(-((x-mu)^2)/(2*sd^2))
      return(p2)
    }
    
    #theta1 is the right starting input for our estimation
    theta1 <- list(c(0,1),c(0,1))
    
    #theta2 is not a list, thus should return an warning
    theta2 <- c(0,1) 
    #theta3 is a list but has not the right amount of parameters, thus should retur an warning
    theta3 <- list(1,c(0,1))
    
    output1 <- multiHMM3(x=z,theta1,m=2,L1,L2)
    

    
    #theta1
    #The estimated thetas should be the same and have a mu =0.5 and delta=0.5 
    expect_equal(output1$Theta[[1]],c(0.5,0.5))
    expect_equal(output1$Theta[[2]],c(0.5,0.5))
    
    
    #compared to the EM-algorithm our code does not produce an error, because the nlminb() function does not find 
    #the right parameters and just returns the input parameters. BUT it does not produce an error. Only a warning.
    
    
    #theta2
    expect_warning(multiHMM3(x=z,theta2,m=2,L1,L2)) 
    #theta3 
    expect_warning(multiHMM3(x=z,theta3,m=2,L1,L2))
})