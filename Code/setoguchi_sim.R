
sim.setoguchi.data <- function(N = 10000, sigma = .1, exposureModel = "E") {
  

## Setoguchi 2008 simulations
    
    
    ## use the original W simulation?
    use.original <- FALSE
    
    ## define the covariance matrix
    Sigma <- diag(1, 10)
    Sigma[1, 5] <- Sigma[5, 1] <- 0.2
    Sigma[2, 6] <- Sigma[6, 2] <- 0.9
    Sigma[3, 8] <- Sigma[8, 3] <- 0.2
    Sigma[4, 9] <- Sigma[9, 4] <- 0.9
    
    
    
    if (!use.original){
      ## simulate the underlying variables
      ## here we use everything continuous
      #set.seed(8794)
      x <- data.frame(mvrnorm(N, mu = rep(0, 10), Sigma))
      colnames(x) <- paste("W.", 1:10, sep = "")
    }
    
    if (use.original){
      ## slightly modified but should still be reasonable
      ## in original simulation 1,3,5,6,8,9 are binary
      W.1 <- 2 * rbinom(N, size = 1, prob = .5)
      W.2 <- rnorm(N)
      W.3 <- 2* rbinom(N, size = 1, prob = .5)
      W.4 <- rnorm(N)
      W.5 <- 2 * rbinom(N, size = 1, prob = .5)
      W.6 <- 2 * rbinom(N, size = 1, prob = .5)
      W.7  <- rnorm(N)
      W.8 <- 2 * rbinom(N, size = 1, prob = .5)
      W.9 <- 2 * rbinom(N, size = 1, prob = .5)
      W.10 <- rnorm(N)
      x <- data.frame(t(chol(Sigma) %*% t(cbind(W.1, W.2 , W.3 , W.4 , W.5 , W.6 , W.7 , W.8 , W.9 , W.10))))
      colnames(x) <- paste("W.", 1:10, sep = "")
    }
    
    ## assign the transformed W's
    W.1 <- x$W.1
    W.2 <- x$W.2
    W.3 <- x$W.3
    W.4 <- x$W.4
    W.5 <- x$W.5
    W.6 <- x$W.6
    W.7 <- x$W.7
    W.8 <- x$W.8
    W.9 <- x$W.9
    W.10 <- x$W.10
    
    ## simulate the exposure model
    if(exposureModel == "A") {
      prob.true <- plogis(
        0 + (0.8 * W.1) - (0.25 * W.2) + (0.6 * W.3) - (0.4 * W.4) - (0.8 * W.5) -
          (0.5 * W.6) + (0.7 * W.7)
      )
    } else if (exposureModel == "E") {
      prob.true <- plogis(
        0 + (0.8 * W.1) - (0.25 * W.2) + (0.6 * W.3) - (0.4 * W.4) - (0.8 * W.5) -
          (0.5 * W.6) + (0.7 * W.7) - (0.25 * W.2 * W.2) + (0.8 * 0.5 * W.1 * W.3) - 
          (0.25 * 0.7 * W.2 * W.4) - (0.4 * 0.5 * W.4 * W.5) - (0.8 * 0.5 * W.5 * W.6)
      )
    } else if (exposureModel == "G") {
      prob.true <- plogis(
        0 + (0.8*W.1) - (0.25*W.2) + (0.6*W.3) - (0.4*W.4) - (0.8*W.5) - (0.5*W.6) + (0.7*W.7) - (0.25*W.2*W.2) -
        (0.4*W.4*W.4) - (0.7*W.7*W.7) + (0.8*0.5*W.1*W.3) + (0.8*0.5*W.1*W.6) - (0.25*0.7*W.2*W.3) - (0.25*0.7*W.2*W.4) +
        (0.6*0.5*W.3*W.4)  + (0.6*0.5*W.3*W.5) - (0.4*0.5*W.4*W.5) - (0.4*0.7*W.4*W.6) - (0.8*0.5*W.5*W.6) - (0.8*0.5*W.5*W.7) 
      )
    } 
    
    
    
  #  trt  <- rbinom(N, size = 1, prob = prob.true)
  #  set.seed(2651)
    x$trt  <- rbinom(N, size = 1, prob = 0.5)
    
    XL2 <- ones <- rep(1,N)
    XL1 <- zeros <- rep(0,N)
 
    XL1[W.3 >= -0.75 & W.3 <=0.75] <- ones[W.3 >= -0.75 & W.3 <=0.75]
    XL2[W.4 >= -0.75 & W.4 <=0.75] <- zeros[W.4 >= -0.75 & W.4 <=0.75]
    
    
    ## simulate the outcome variable
    ## exposureModels E and G are modified using Stuart PSMG December 2015
    eps <- rnorm(N, sd = sigma)
    
    
    if(exposureModel == "A") {
      x$y0.true <- -(-3.85 + 0.3*W.1 - 0.36*W.2 - 0.73*W.3 - 0.2*W.4 + 0.71*W.8 - 0.19*W.9 +
                     0.26*W.10 - 0.4*0)
      
      x$y1.true <- -(-3.85 + 0.3*W.1 - 0.36*W.2 - 0.73*W.3 - 0.2*W.4 + 0.71*W.8 - 0.19*W.9 +
                     0.26*W.10 - 0.4*1)
      
      x$y <- -(-3.85 + 0.3*W.1 - 0.36*W.2 - 0.73*W.3 - 0.2*W.4 + 0.71*W.8 - 0.19*W.9 +
               0.26*W.10 - 0.4*x$trt + eps)
      
      
    } else if (exposureModel == "E") {
      
      x$y0.true <- -(-3.85 + 0.3*W.1 - 0.36*W.2 - 0.73*W.3 - 0.2*W.4 + 0.71*W.8 - 0.19*W.9 +
                     0.26*W.10 - 0.4*0 - 0.36*0.5*W.2*W.2 - 0.73*0.7*W.1*W.3 - 0.2*W.1*W.4)
      
      x$y1.true <- -(-3.85 + 0.3*W.1 - 0.36*W.2 - 0.73*W.3 - 0.2*W.4 + 0.71*W.8 - 0.19*W.9 +
                     0.26*W.10 - 0.4*1 - 0.36*0.5*W.2*W.2 - 0.73*0.7*W.1*W.3 - 0.2*W.1*W.4)
      
      x$y <- -(-3.85 + 0.3*W.1 - 0.36*W.2 - 0.73*W.3 - 0.2*W.4 + 0.71*W.8 - 0.19*W.9 +
               0.26*W.10 - 0.4*x$trt - 0.36*0.5*W.2*W.2 - 0.73*0.7*W.1*W.3 - 0.2*W.1*W.4 + eps)
      
      
      
    } else if (exposureModel == "G") {
      
      x$y0.true <- -(-3.85 +0.3*W.1 - 0.36*W.2 - 0.73*W.3 - 0.2*W.4 + 0.71*W.8 - 0.19*W.9 +
                     0.26*W.10 - 0.4*0 - 0.19*W.7*W.6 - 0.19*W.9*W.2*W.3 + 0.26*W.10*W.1 +
                     0.71*0.7*W.8*W.4 - 0.19*0.3*W.9*W.2*W.3 + 0.26*W.10*W.1 - 0.36*0.5*W.2*W.2 -
                     0.73*W.3*W.3 - 0.19*W.6*W.6 + .8*W.1*W.5*XL1*0 - .4*W.2*W.6*XL2*0 - .2*W.7*W.8*XL1)
      
      x$y1.true <- -(-3.85 +0.3*W.1 - 0.36*W.2 - 0.73*W.3 - 0.2*W.4 + 0.71*W.8 - 0.19*W.9 +
                     0.26*W.10 - 0.4*1 - 0.19*W.7*W.6 - 0.19*W.9*W.2*W.3 + 0.26*W.10*W.1 +
                     0.71*0.7*W.8*W.4 - 0.19*0.3*W.9*W.2*W.3 + 0.26*W.10*W.1 - 0.36*0.5*W.2*W.2 -
                     0.73*W.3*W.3 - 0.19*W.6*W.6 + .8*W.1*W.5*XL1*1 - .4*W.2*W.6*XL2*1 - .2*W.7*W.8*XL1)
      
      x$y      <-  -(-3.85 +0.3*W.1 - 0.36*W.2 - 0.73*W.3 - 0.2*W.4 + 0.71*W.8 - 0.19*W.9 +
                     0.26*W.10 - 0.4*x$trt - 0.19*W.7*W.6 - 0.19*W.9*W.2*W.3 + 0.26*W.10*W.1 +
                     0.71*0.7*W.8*W.4 - 0.19*0.3*W.9*W.2*W.3 + 0.26*W.10*W.1 - 0.36*0.5*W.2*W.2 -
                     0.73*W.3*W.3 - 0.19*W.6*W.6 + .8*W.1*W.5*XL1*x$trt - .4*W.2*W.6*XL2*x$trt - .2*W.7*W.8*XL1  + 
                     eps)
    } 
    

 
  x$TE = x$y1.true - x$y0.true
  
  ## assemble the data as a bivariate missing data problem
  x$y0 <- x$y1 <- rep(NA, N)
  x$y0[x$trt == 0] <- x$y[x$trt == 0]
  x$y1[x$trt == 1] <- x$y[x$trt == 1]
  
  ## return the simulated data
  data.frame(x)
  
}