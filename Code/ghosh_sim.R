
#Ghosh Simulation

sim.ghosh.data <- function(N = 10000, sigma = .1, case = 1, linearF = TRUE) {
  
  
  ## simulate x and treatment variable
  set.seed(111)
  x.c <- matrix(rnorm(N * 11), N, 11)
  
  set.seed(222)
  x.d <- matrix(rbinom(N * 9, size = 1, prob = .5), N, 9)
  
  x <- data.frame(x.c, x.d)
  colnames(x) <- paste("X",1:ncol(x),sep="")
  
  if (linearF) {##simulation extended to allow nonlinear F
    F <- (-2.0 + .028 * x[, 1] - 0.374 * x[, 2] -.03 * x[, 3] + .118 * x[, 4] 
          - 0.394 * x[, 11] + .875 * x[, 12] + .9 * x[, 13])
  }
  else {
    ## nonlinear F
    F <- (-2.0 + .028 * x[, 1] - 3.374 * x[, 2] -.03 * x[, 3] + .118 * x[, 4] ^ 2 
          - 1.394 * x[, 11] + 1.875 * x[, 7] * x[, 8])
  }
  prob.true <- plogis(F)
  
  set.seed(454)
  x$trt <- rbinom(N, size = 1, prob = prob.true)
  
  ## noise
  eps <- rnorm(N, sd = sigma)
  
  ## simulate the outcome
  if (case == 1) {##ghosh simulation case=1: no confounding
    x$y.cntrl.true <- (2.455 + .40 * 0 + .1 * x[, 1] - .154 * x[, 2]
                       - .152 * x[, 11] - .126 * x[, 12])
    x$y.trt.true <- (2.455 + .40 * 1 + .1 * x[, 1] - .154 * x[, 2]
                     - .152 * x[, 11] - .126 * x[, 12])
    x$y <- (2.455 + .40 * x$trt + .1 * x[, 1] - .154 * x[, 2]
            - .152 * x[, 11] - .126 * x[, 12] + eps)    
    
  }
  else if (case == 2) {##ghosh simulation case=2: X11 is confounded variable
    x$y.cntrl.true <- (2.455 + .40 * 0 + .1 * x[, 1] - .154 * x[, 2]
                       - .152 * x[, 11] - .126 * x[, 12] - .3 * 0 * x[, 11])
    x$y.trt.true <- (2.455 + .40 * 1 + .1 * x[, 1] - .154 * x[, 2]
                     - .152 * x[, 11] - .126 * x[, 12] - .3 * 1 * x[, 11])
    x$y          <-  (2.455 + .40 * x$trt + .1 * x[, 1] - .154 * x[, 2]
                      - .152 * x[, 11] - .126 * x[, 12] - .3 * x$trt * x[, 11]
                      + eps)
  }
  else if (case == 3) {#ghosh case=2 modified with additional confounded variable (X2)
    x$y.cntrl.true <- (2.455 + .40 * 0 + .1 * x[, 1] - .154 * x[, 2]
                       - .152 * x[, 11] - .126 * x[, 12] - .3 * 0 * x[, 11])
    x$y.trt.true <- (2.455 + .40 * 1 + .1 * x[, 1] - .254 * x[, 2] ^ 2
                     - .152 * x[, 11] - .126 * x[, 12] - .3 * 1 * x[, 11])
    x$y <- (2.455 + .40 * x$trt + .1 * x[, 1]
            - .154 * x[, 2] * (1 - x$trt) - .254 * x[, 2] ^ 2 * x$trt
            - .152 * x[, 11] - .126 * x[, 12] - .3 * x$trt * x[, 11]
            + eps)
    
  }
  
  x$TE = x$y.trt.true - x$y.cntrl.true
  
  x$y1 <- x$y0 <- rep(NA, N)
  x$y0[x$trt == 0] <- x$y[x$trt == 0]
  x$y1[x$trt == 1] <- x$y[x$trt == 1]
  
# if we dont want noise, then instead of taking control & trt from y, we take them from y.cntrl.true and y.trt.true  
#  x$y.cntrl.true = ifelse(x$trt==0 , x$y.cntrl.true , NA)
#  x$y.trt.true = ifelse(x$trt==1 , x$y.trt.true , NA)

    return(x)
}