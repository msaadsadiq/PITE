
#importing the simulation functions
source("~/Dropbox/BioStats/PITE/andrea_sim.R")
source("~/Dropbox/BioStats/PITE/ghosh_sim.R")
source("~/Dropbox/BioStats/PITE/setoguchi_sim.R")


#Variables

nrep <- c(250)
sim.type <- c("Andrea", "Ghosh","Setoguchi")[3]
forest.type <- c("randomForest","rfsrc","rfsrcsyn")[3]

## forest parameters
ntree            <- c(10, 50, 100, 250, 500, 1000)[6]
nsplit           <- c(0, 1, 2, 10)[3]
nodesize.grw     <- c(1, 3, 5, 10)[2]

## synthetic forests
nodesizeSeq      <- c(1:10, 20, 30, 50, 100)
mtrySeq          <- c(1, 10, 20)

## verbose flag
verbose <- FALSE 


#monte-carlo simulation
rO <- lapply(1:nrep, function(b) {
  
  ## Display current loop index
  if (nrep > 1) {
    cat("Replication:", b, "\n")
  }
  
  if(sim.type == "Andrea") {
    data <- sim.andrea.data(10000,7)
    
    estimator.cntrl.data <- data.frame(data[1:5000,1:7])
    estimator.cntrl.data$control <-data[1:5000,78]
    
    test.cntrl.data <- data.frame(data[5001:10000,1:7])
    test.cntrl.data$control <-data[5001:10000,78]
    
    estimator.trt.data <- data.frame(data[1:5000,1:7])
    estimator.trt.data$trt <-data[1:5000,79]
    
    test.trt.data <- data.frame(data[5001:10000,1:7])
    test.trt.data$trt <-data[5001:10000,79]
    
  } else if(sim.type == "Ghosh") {
    data <- sim.ghosh.data(N = 10000, sigma = .1, linearF = TRUE)
    
    estimator.cntrl.data <- data.frame(data[1:5000,1:20])
    estimator.cntrl.data$control <-data[1:5000,26]
    
    test.cntrl.data <- data.frame(data[5001:10000,1:20])
    test.cntrl.data$control <-data[5001:10000,26]
    
    estimator.trt.data <- data.frame(data[1:5000,1:20])
    estimator.trt.data$trt <-data[1:5000,27]
    
    test.trt.data <- data.frame(data[5001:10000,1:20])
    test.trt.data$trt <-data[5001:10000,27]
    
  } else if(sim.type == "Setoguchi") {
    data <- sim.setoguchi.data(N = 10000, sigma = .1, exposureModel = "G")
    
    estimator.cntrl.data <- data.frame(data[1:5000,1:10])
    estimator.cntrl.data$control <-data[1:5000,17]
    
    test.cntrl.data <- data.frame(data[5001:10000,1:10])
    test.cntrl.data$control <-data[5001:10000,17]
    
    estimator.trt.data <- data.frame(data[1:5000,1:10])
    estimator.trt.data$trt <-data[1:5000,16]
    
    test.trt.data <- data.frame(data[5001:10000,1:10])
    test.trt.data$trt <-data[5001:10000,16]
    
    
  }

########################################################
# First predictive model using Multiple Imputations
########################################################

# Transform data before imputing
  mdf.cntrl <- missing_data.frame(estimator.cntrl.data)
  mdf.cntrl <- change(mdf.cntrl, y= "control", what ="transformation", to="identity")
#  mdf.cntrl <- change(mdf.cntrl, y= "control", what ="imputation_method", to="mean")  
  mdf.trt <- missing_data.frame(estimator.trt.data)
  mdf.trt <- change(mdf.trt, y= "trt", what ="transformation", to="identity")
#  mdf.cntrl <- change(mdf.trt, y= "trt", what ="imputation_method", to="mean")  
  
# Impute missing values using mi in control column
 
  mi.cntrl.impute  <- mi(mdf.cntrl, n.iter=100, n.chains=100, save_models = TRUE)
  mi.cntrl.pred <- sapply(complete(mi.cntrl.impute, m=100), function(x) {x$control})
  
# Impute missing values by mi in treatment column
 
  mi.trt.impute  <- mi(mdf.trt, n.iter=100, n.chains =100, save_models = TRUE)
  mi.trt.pred <- sapply(complete(mi.trt.impute, m=100), function(x) {x$trt})
 
# Evaluating Coefficients for both trt & cntrl
  
  if(sim.type =="Andrea") {
    pool.cntrl.coefficients <-  pool(control ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, mi.cntrl.impute, m=100)
    pool.trt.coefficients <- pool(trt ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, mi.trt.impute, m=100)
  } else if(sim.type =="Ghosh") {
    pool.cntrl.coefficients <-  pool(control ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20, mi.cntrl.impute, m=100)
    pool.trt.coefficients <- pool(trt ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, mi.trt.impute, m=100)
  } else if(sim.type =="Setoguchi") {
    pool.cntrl.coefficients <-  pool(control ~ W.1 + W.2 + W.3 + W.4 + W.5 + W.6 + W.7 + W.8 + W.9 + W.10, mi.cntrl.impute, m=100)
    pool.trt.coefficients   <-  pool(trt ~ W.1 + W.2 + W.3 + W.4 + W.5 + W.6 + W.7 + W.8 + W.9 + W.10, mi.trt.impute, m=100)
  }
  
  cntrl.coefficients <- summary(pool.cntrl.coefficients)$coefficients[2:11,1]
  trt.coefficients   <- summary(pool.trt.coefficients)$coefficients[2:11,1]
  
# PITE and Bias for Multiple Imputation
  train.mi.pite <- rowMeans(mi.trt.pred - mi.cntrl.pred)
  train.mi.bias <- train.mi.pite - data$TE[1:5000]
  train.mi.bias.rel <- train.mi.bias / data$TE[1:5000]
  
  
  
#######################################################
# Second predictive model using Random Forests
#######################################################
  
if(forest.type == "randomForest") {

  #using random forests
    rf.cntrl.est <- randomForest(control ~ ., estimator.cntrl.data,  nodesize=100, na.action=na.roughfix)
    estimator.cntrl.data$pred.cntrl <- rf.cntrl.est$predicted
  
  #2 find missing values by mi in treatment column
    rf.trt.est <- randomForest(trt ~ ., estimator.trt.data,  nodesize=12, na.action=na.roughfix)
    estimator.trt.data$pred.trt <- rf.trt.est$predicted
  
} else if(forest.type == "rfsrc") {
 
  # Regress control column values 
    rf.cntrl.est <- rfsrc(control ~ ., estimator.cntrl.data,  importance = "none", 
                        ntree = ntree, nsplit = nsplit, nodesize = nodesize.grw)
  
  # Regress trt column values
    rf.trt.est <- rfsrc(trt ~ ., estimator.trt.data,  importance = "none", 
                      ntree = ntree, nsplit = nsplit, nodesize = nodesize.grw)

} else if(forest.type == "rfsrcsyn") {

  # Regress control column values     
  rf.cntrl.est  <-   rfsrcSyn(control ~ ., estimator.cntrl.data, importance = "none",
                                            ntree = ntree, nsplit = nsplit, nodesizeSeq = nodesizeSeq,
                                             mtrySeq = mtrySeq, nodesize = nodesize.grw, 
                                            verbose = verbose)

  # Regress trt column values     
  rf.trt.est    <-   rfsrcSyn(trt ~ ., estimator.trt.data, importance = "none",
                                            ntree = ntree, nsplit = nsplit, nodesizeSeq = nodesizeSeq, 
                                            mtrySeq = mtrySeq, nodesize = nodesize.grw, 
                                            verbose = verbose)
}
  
  
  estimator.cntrl.data$pred.rf     <-rfsrcSyn(object = rf.trt.est, newdata = estimator.cntrl.data[,1:10])$rfSynPred$predicted 
  estimator.cntrl.data$pred.rf.opp <-rfsrcSyn(object = rf.cntrl.est, newdata = estimator.cntrl.data[,1:10])$rfSynPred$predicted
  estimator.trt.data$pred.rf       <-rfsrcSyn(object = rf.cntrl.est, newdata = estimator.trt.data[,1:10])$rfSynPred$predicted
  estimator.trt.data$pred.rf.opp   <-rfsrcSyn(object = rf.trt.est, newdata = estimator.trt.data[,1:10])$rfSynPred$predicted    
  
  # Test data PITE RF & Bias
  train.rf.pite <- estimator.trt.data$pred.rf - estimator.cntrl.data$pred.rf
  train.rf.bias <- train.rf.pite - data$TE[1:5000]
  train.rf.bias.rel <- train.rf.bias / data$TE[1:5000]
  
  # opposite forests
  train.rf.pite.opp <-estimator.trt.data$pred.rf.opp - estimator.cntrl.data$pred.rf.opp
  train.rf.bias.opp <- train.rf.pite.opp - data$TE[1:5000]
  train.rf.bias.rel.opp <- train.rf.bias.opp / data$TE[1:5000]
  
  # RMSE
  data.train <- data$y1.true[1:5000] - data$y0.true[1:5000]
  train.RMSE.rf <- sqrt( (train.rf.pite - data.train) ^2 )
  train.RMSE.mi <- sqrt( (train.mi.pite - data.train) ^2 )
  
  
  
#######################################################
# Using Models to predict testing data
####################################################### 
 
########  MI  ######### 
 
 test.cntrl.data$pred.mi <-  cntrl.coefficients[1] * test.cntrl.data$W.1 + 
                                   cntrl.coefficients[2] * test.cntrl.data$W.2 + 
                                   cntrl.coefficients[3] * test.cntrl.data$W.3 + 
                                   cntrl.coefficients[4] * test.cntrl.data$W.4 + 
                                   cntrl.coefficients[5] * test.cntrl.data$W.5 + 
                                   cntrl.coefficients[6] * test.cntrl.data$W.6 + 
                                   cntrl.coefficients[7] * test.cntrl.data$W.7 +  
                                   cntrl.coefficients[8] * test.cntrl.data$W.8 + 
                                   cntrl.coefficients[9] * test.cntrl.data$W.9 + 
                                   cntrl.coefficients[10] * test.cntrl.data$W.10 
 
 test.trt.data$pred.mi <-  trt.coefficients[1] * test.trt.data$W.1 + 
                               trt.coefficients[2] * test.trt.data$W.2 + 
                               trt.coefficients[3] * test.trt.data$W.3 + 
                               trt.coefficients[4] * test.trt.data$W.4 + 
                               trt.coefficients[5] * test.trt.data$W.5 + 
                               trt.coefficients[6] * test.trt.data$W.6 + 
                               trt.coefficients[7] * test.trt.data$W.7 + 
                               trt.coefficients[8] * test.trt.data$W.8 + 
                               trt.coefficients[9] * test.trt.data$W.9 + 
                               trt.coefficients[10] * test.trt.data$W.10
   
 # Test data PITE MI 
 test.mi.pite <- test.trt.data$pred.mi - test.cntrl.data$pred.mi
 test.mi.bias <- test.mi.pite - data$TE[5001:10000]
 test.mi.bias.rel <- test.mi.bias / data$TE[1:5000]
 
 ########  RF  ######### 
 
 # predicting the model
 
 if(forest.type == "rfsrc") {
   test.cntrl.data$pred.rf <- predict(rf.cntrl.est, test.cntrl.data[,1:10])$predicted
   test.trt.data$pred.rf <- predict(rf.trt.est, test.trt.data[,1:10])$predicted 
  
 } 
 
 else if(forest.type == "rfsrcsyn") {
   test.cntrl.data$pred.rf <- rfsrcSyn(object = rf.cntrl.est, newdata = test.cntrl.data[,1:10])$rfSynPred$predicted
   test.trt.data$pred.rf   <- rfsrcSyn(object = rf.trt.est, newdata = test.trt.data[,1:10])$rfSynPred$predicted
 }
 
 # Test data PITE RF & Bias
 test.rf.pite <- test.trt.data$pred.rf - test.cntrl.data$pred.rf
 test.rf.bias <- test.rf.pite - data$TE[5001:10000]
 test.rf.bias.rel <- test.rf.bias / data$TE[1:5000]
 
 
 ########  Variability ######
 
 variance.mi <- (test.mi.pite - mean(test.mi.pite))^2
 variance.rf <- (test.rf.pite - mean(test.rf.pite))^2
 
 ########  RMSE  ######### 
 data.test <- data$y1.true[5001:10000] - data$y0.true[5001:10000]
 test.RMSE.rf <- sqrt( (test.rf.pite - data.test) ^2 )
 test.RMSE.mi <- sqrt( (test.mi.pite - data.test) ^2 )
 
 
#### return the goodies
list(
  train.rf.pite = train.rf.pite,           #A
  train.mi.pite = train.mi.pite,           #B
  test.rf.pite  = test.rf.pite,            #C
  test.mi.pite  = test.mi.pite,            #D
  train.rf.bias = train.rf.bias,           #E
  train.mi.bias = train.mi.bias,           #F
  test.rf.bias  = test.rf.bias,            #G
  test.mi.bias  = test.mi.bias,            #H
  TE            = data$TE,                 #I
  train.RMSE.rf = train.RMSE.rf,           #J
  train.RMSE.mi = train.RMSE.mi,           #K
  test.RMSE.rf  = test.RMSE.rf,            #L
  test.RMSE.mi  = test.RMSE.mi,            #M
  variance.rf   = variance.rf,             #N
  variance.mi   = variance.mi,             #O
  train.rf.bias.rel = train.rf.bias.rel,   #P
  train.mi.bias.rel = train.mi.bias.rel,   #Q
  test.rf.bias.rel = test.rf.bias.rel,     #R
  test.mi.bias.rel = test.mi.bias.rel,      #S
  train.rf.pite.opp = train.rf.pite.opp,    #T
  train.rf.bias.opp = train.rf.bias.opp,    #U
  train.rf.bias.rel.opp = train.rf.bias.rel.opp #V 
  
   ) 

})

#### Label list of lists for easy reference
rO <- lapply(rO,FUN=function(x) {names(x) <-   
  c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V")
x 
})




### Average each component of the monte-carlo simulation

avg.train.rf.pite  <-  rowMeans(sapply(rO, function(x) {x$A}), na.rm = TRUE)
avg.train.mi.pite  <-  rowMeans(sapply(rO, function(x) {x$B}), na.rm = TRUE)

avg.test.rf.pite   <-  rowMeans(sapply(rO, function(x) {x$C}), na.rm = TRUE)
avg.test.mi.pite   <-  rowMeans(sapply(rO, function(x) {x$D}), na.rm = TRUE)

avg.train.rf.bias  <-  rowMeans(sapply(rO, function(x) {x$E}), na.rm = TRUE)
avg.train.mi.bias  <-  rowMeans(sapply(rO, function(x) {x$F}), na.rm = TRUE)

avg.test.rf.bias   <-  rowMeans(sapply(rO, function(x) {x$G}), na.rm = TRUE)
avg.test.mi.bias   <-  rowMeans(sapply(rO, function(x) {x$H}), na.rm = TRUE)

meanTE             <- rowMeans(sapply(rO, function(x) {x$I}), na.rm = TRUE)

avg.train.rmse.rf <- rowMeans(sapply(rO, function(x) {x$J}), na.rm = TRUE)
avg.train.rmse.mi <- rowMeans(sapply(rO, function(x) {x$K}), na.rm = TRUE)

avg.test.rmse.rf <- rowMeans(sapply(rO, function(x) {x$L}), na.rm = TRUE)
avg.test.rmse.mi <- rowMeans(sapply(rO, function(x) {x$M}), na.rm = TRUE)

avg.variability.rf <-  rowMeans(sapply(rO, function(x) {x$N}), na.rm = TRUE)
avg.variability.mi <-  rowMeans(sapply(rO, function(x) {x$O}), na.rm = TRUE)

avg.train.rf.bias.rel <-  rowMeans(sapply(rO, function(x) {x$P}), na.rm = TRUE)
avg.train.mi.bias.rel <-  rowMeans(sapply(rO, function(x) {x$Q}), na.rm = TRUE)

avg.test.rf.bias.rel <-  rowMeans(sapply(rO, function(x) {x$R}), na.rm = TRUE)
avg.test.mi.bias.rel <-  rowMeans(sapply(rO, function(x) {x$S}), na.rm = TRUE)

avg.train.rf.pite.opp <-  rowMeans(sapply(rO, function(x) {x$T}), na.rm = TRUE)
avg.train.rf.bias.opp <-  rowMeans(sapply(rO, function(x) {x$U}), na.rm = TRUE)
avg.train.rf.bias.rel.opp <-  rowMeans(sapply(rO, function(x) {x$V}), na.rm = TRUE)


#######################################################
# PLOTS & Results
#######################################################  


### Individual Level Bias comparison between MI, RF
par(mfrow=c(1,2))

plot(density(avg.train.rf.bias), xlim=c(-10,10), ylim=c(0,1.5), xlab="Mean Bias", col="red", main="Distribution of individual-level bias - Training")
lines(density(avg.train.mi.bias), col="blue")
lines(density(avg.train.rf.bias.opp), col="green")
legend(-10,1,c("RF","MI","Opp"),lty=c(1,1),lwd=c(1.5,2.5),col=c("red","blue","green"))

plot(density(avg.test.rf.bias), xlim=c(-10,10), ylim=c(0,1.5),xlab="Mean Bias", col="red", main="Distribution of individual-level bias - Testing")
lines(density(avg.test.mi.bias), col="blue")
legend(-10,1,c("RF","MI"),lty=c(1,1),lwd=c(1.5,2.5),col=c("red","blue"))


### Individual Level Relative Bias comparison between MI, RF
par(mfrow=c(1,2))

plot(density(avg.train.rf.bias.rel),  xlim=c(-300,100), ylim=c(0,2),xlab="Mean Bias", col="red", main="Distribution of Rel. bias - Training")
lines(density(avg.train.mi.bias.rel), col="blue")
lines(density(avg.train.rf.bias.rel.opp), col="green")
legend(-200,1,c("RF","MI","Opp"),lty=c(1,1),lwd=c(1.5,2.5),col=c("red","blue","green"))

plot(density(avg.test.rf.bias.rel),  xlim=c(-300,100), ylim=c(0,2), xlab="Mean Bias", col="red", main="Distribution of Rel. bias - Testing")
lines(density(avg.test.mi.bias.rel), col="blue")
legend(-200,1,c("RF","MI"),lty=c(1,1),lwd=c(1.5,2.5),col=c("red","blue"))



### Bias as a function of True Treatment Effect
par(mfrow=c(2,2))
plot(meanTE[1:5000], avg.train.mi.bias,xlim=c(-7,7), ylim=c(-15,10), xlab="True treatment effect", ylab="Bias - MI", main="mi Bias vs true treatment effects - Training")
plot(meanTE[1:5000], avg.train.rf.bias,xlim=c(-7,7), ylim=c(-15,10), xlab="True treatment effect", ylab="Bias - RF",main="rf Bias vs true treatment effects - Training")

plot(meanTE[5001:10000], avg.test.mi.bias,xlim=c(-7,7), ylim=c(-15,10), xlab="True treatment effect", ylab="Bias - MI", main="mi Bias vs true treatment effects - Testing")
plot(meanTE[5001:10000], avg.test.rf.bias,xlim=c(-7,7), ylim=c(-15,10), xlab="True treatment effect", ylab="Bias - RF",main="rf Bias vs true treatment effects - Testing")

# Comparison of True and Predicted values under treatment and control 
par(mfrow=c(2,2))

plot(estimator.cntrl.data$pred.rf, data$y0.true[1:5000], xlab="Predicted Y - rfsrc", ylab="True Y - control condition")
plot(estimator.trt.data$pred.rf  , data$y1.true[1:5000], xlab="Predicted Y - rfsrc", ylab="True Y - treatment condition")
plot(rowMeans(mi.cntrl.pred),      data$y0.true[1:5000], xlab="Predicted Y - mi", ylab="True Y - control condition")
plot(rowMeans(mi.trt.pred)  ,      data$y1.true[1:5000], xlab="Predicted Y - mi", ylab="True Y - treatment condition")

plot(test.cntrl.data$pred.rf, data$y0.true[5001:10000], xlab="Predicted Y - rfsrc", ylab="True Y - control condition")
plot(test.trt.data$pred.rf, data$y1.true[5001:10000], xlab="Predicted Y - rfsrc", ylab="True Y - treatment condition")
plot(test.cntrl.data$pred.mi, data$y0.true[5001:10000], xlab="Predicted Y - mi", ylab="True Y - control condition")
plot(test.trt.data$pred.mi, data$y1.true[5001:10000], xlab="Predicted Y - mi", ylab="True Y - treatment condition")



#  True vs Predicted Treatment Effect 
par(mfrow=c(2,2))

plot(avg.train.mi.pite, meanTE[1:5000], xlim=c(-10,10), ylim=c(-6,6), xlab="PITE MI",ylab="True Effect", main="True vs Predicted - Training")
plot(avg.train.rf.pite, meanTE[1:5000], xlim=c(-10,10), ylim=c(-6,6),xlab="PITE RF",ylab="True Effect", main="True vs Predicted - Training")

plot(avg.test.mi.pite, meanTE[5001:10000], xlim=c(-10,10), ylim=c(-6,6),xlab="PITE MI",ylab="True Effect", main="True vs Predicted - Testing")
plot(avg.test.rf.pite, meanTE[5001:10000],xlim=c(-10,10), ylim=c(-6,6), xlab="PITE RF",ylab="True Effect", main="True vs Predicted - Testing")


# Variance of MI and RF estimators vs True treatment effect
par(mfrow=c(1,2))

plot(meanTE[5001:10000], avg.variability.mi, ylim=c(0,10), xlab="True Treatment Effect", ylab="Variance MI")
plot(meanTE[5001:10000], avg.variability.rf, ylim=c(0,10), xlab="True Treatment Effect", ylab="Variance RF")


# Distribution of the RMSE across individuals
par(mfrow=c(1,2))

plot(density(avg.train.rmse.rf),xlim=c(-1,10), ylim=c(0,2), xlab="RMSE",ylab="Density", col="red", main="Dist of RMSE - Training")
lines(density(avg.train.rmse.mi), col="blue")
legend(4,1.2,c("RF","MI"),lty=c(1,1),lwd=c(2.5,2.5),col=c("red","blue"))

plot(density(avg.test.rmse.rf),xlim=c(-1,10), ylim=c(0,2),xlab="RMSE",ylab="Density", col="red", main="Dist of RMSE - Testing")
lines(density(avg.test.rmse.mi), col="blue")
legend(3,1.5,c("RF","MI"),lty=c(1,1),lwd=c(2.5,2.5),col=c("red","blue"))


#  RMSE vs PITE
par(mfrow=c(2,2))

plot(avg.train.rf.pite, avg.train.rmse.rf, xlim=c(-10,10), ylim=c(0,12), xlab="PITE-RF",ylab="RMSE RF",  main=" RMSE vs PITE - Training")
plot(avg.train.mi.pite, avg.train.rmse.mi, xlim=c(-10,10), ylim=c(0,12), xlab="PITE-MI",ylab="RMSE MI",  main=" RMSE vs PITE - Training")
plot(avg.test.rf.pite, avg.test.rmse.rf,  xlim=c(-10,10), ylim=c(0,12), xlab="PITE-RF",ylab="RMSE RF",    main=" RMSE vs PITE - Testing")
plot(avg.test.mi.pite, avg.test.rmse.mi,  xlim=c(-10,10), ylim=c(0,12), xlab="PITE-MI",ylab="RMSE MI",    main=" RMSE vs PITE - Testing")


#  Rel. Bias vs PITE
par(mfrow=c(2,2))

plot(avg.train.rf.pite, avg.train.rf.bias.rel, xlim=c(-8,10), ylim=c(0,50), xlab="PITE-RF",ylab="Rel.Bias RF",  main=" Rel.Bias vs PITE - Training")
plot(avg.train.mi.pite, avg.train.mi.bias.rel, xlim=c(-8,10), ylim=c(0,50), xlab="PITE-MI",ylab="Rel.Bias MI",  main=" Rel.Bias vs PITE - Training")
plot(avg.test.rf.pite, avg.test.rf.bias.rel,  xlim=c(-8,10), ylim=c(0,50), xlab="PITE-RF",ylab="Rel.Bias RF",    main=" Rel.Bias vs PITE - Testing")
plot(avg.test.mi.pite, avg.test.mi.bias.rel,  xlim=c(-8,10), ylim=c(0,50), xlab="PITE-MI",ylab="Rel.Bias MI",    main=" Rel.Bias vs PITE - Testing")



#  Bias vs PITE
par(mfrow=c(2,2))

plot(avg.train.rf.pite, avg.train.rf.bias, xlim=c(-5,10), ylim=c(0,10), xlab="PITE-RF",ylab="Bias RF",  main=" Bias vs PITE - Training")
plot(avg.train.mi.pite, avg.train.mi.bias, xlim=c(-5,10), ylim=c(0,10), xlab="PITE-MI",ylab="Bias MI",  main=" Bias vs PITE - Training")
plot(avg.test.rf.pite, avg.test.rf.bias,  xlim=c(-5,10), ylim=c(0,10), xlab="PITE-RF",ylab="Bias RF",    main="Bias vs PITE - Testing")
plot(avg.test.mi.pite, avg.test.mi.bias,  xlim=c(-5,10), ylim=c(0,10), xlab="PITE-MI",ylab="Bias MI",    main="Bias vs PITE - Testing")


















































  
 
  
