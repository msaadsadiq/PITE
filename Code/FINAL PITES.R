#LOAD LIBRARIES#

library(ggplot2)
library(snow) #for parallelization
library(mi)

#FUNCTION TO GENERATE ONE SIM#
dataGen=function(N,p,seed){
 set.seed(2398)
 X=rbinom(N*p,1,.5)
 df=data.frame(matrix(X,nrow=N,ncol=p))
 df$obs.txt=rep(0:1,N/2)
 df$TE=-1.3*df$X1-1.2*df$X2-.6*df$X3+.3*df$X4+.5*df$X5+1.1*df$X6+1.2*df$X7

 seed=set.seed(seed)
 df$y.obs.control=rnorm(N,0,1)#observed y value under control
 df$y.obs.tx= ifelse(df$obs.txt==1, (df$y.obs.control + df$TE),NA) #observed y value under TX
 #df$Y=ifelse(df$obs.txt==0,df$y.obs.control,df$y.obs.tx) #observed Y value
 df$y.obs.control=ifelse(df$obs.txt==0,df$y.obs.control,NA) #observed y value under control
 df$ob=rep(0:1,each=N/2)
 df$sim=rep(length(seed),each=N)
 return(df)}

#CREATE MI FUNCTION#
mi.func<-function(x){
txt.imp=as.data.frame(x[,c(79,1:75)])
cont.imp=as.data.frame(x[,c(78,1:75)])
names(cont.imp)
#Confirms the imputation formula
info.i<-mi.info(txt.imp) 
info.i[["y.obs.tx"]]$imp.formula
info.c.i<-mi.info(cont.imp) 
info.c.i[["y.obs.control"]]$imp.formula

#This does the imputation for the control condition (imputes the treatment)
mi.control.i<-mi (txt.imp, info.i, n.imp=nimp, max.minutes=20000)
midf.cont.i=mi.completed(mi.control.i)

conv.c=mi.control.i@converged
c.m=(mi.control.i@m)
conv.m.c=cbind(conv.c, c.m)

#This does the imputation for the treatment condition (imputes the control)
mi.txt.i<-mi (cont.imp, info.c.i, n.imp=nimp, max.minutes=20000)
midf.txt.i=mi.completed(mi.txt.i)

conv.t=mi.txt.i@converged
t.m=(mi.txt.i@m)
conv.m.t=cbind(conv.t, t.m)


#Get coefficients for each chain (imputation)
tf <- vector("list", nimp)####change number of imputes
for (j in 1:nimp){ #########CHANGE
        s <- paste("Chain",j, sep="")
        tf[[j]] <- mi.control.i@imp[[s]]$y.obs.tx@model$coefficients
}
tf
c=do.call(rbind,tf) #converts list to matrix

cf <- vector("list", nimp)####change number of imputes
for (z in 1:nimp){ #############CHANGE
        s <- paste("Chain",z, sep="")
        cf[[z]] <- mi.txt.i@imp[[s]]$y.obs.control@model$coefficients
}
cf #Each row is a set of coefficients for that imputation
b=do.call(rbind,cf)

parms=list(c,b,conv.m.c,conv.m.t)
return(parms)
}


system.time()

#RUN ANALYSES
Nsims=100
  #number of simulatins
N=10000
  #number of cases per sim
p=75
  #number of predictors 
nimp=100
  #number of imputations per sim
predictors=c(1:75)
  #predictors that go into the model
set.seed(740)
seedlist=round(runif(Nsims,0,10000))

## Step 1: Create a cluster of child processes 
cl <- makeCluster( 10 )

## Step 2: Load the necessary R package(s)
## N.B. length(cl) is the number of child processes
##      in the cluster 

par.setup <- parLapply( cl, 1:length(cl),
    function(x) {
     	require(mi)})

## Step 3: Distribute the necessary R objects 
clusterExport( cl, c('seedlist','dataGen', 'mi.func','Nsims','N','p','nimp','predictors'))

## Step 4: Do the computation
par.output <- parLapply(cl, seedlist,
    function(x) {
       filename=x
 data=dataGen(N,p,x) #generate 100 cases with 75 predictors with seeds = x (aka "seedlist")
mi= mi.func(data)
coeff.t=(matrix(unlist(mi[1]),nrow=76,byrow=T))
coeff.c=(matrix(unlist(mi[2]),nrow=76,byrow=T))
    

#Raw data
raw.ob=as.matrix(data[which(data$ob==1),1:75])
raw.ob=as.matrix(cbind(rep(1,(dim(raw.ob)[1])),raw.ob))
raw.in=as.matrix(data[which(data$ob==0),1:75])
raw.in=as.matrix(cbind(rep(1,(dim(raw.in)[1])),raw.in))

#Calculating predicted values
	#Control
pred.c.ob=raw.ob%*%(coeff.c) #outbag
pred.c.in=raw.in%*%(coeff.c)#insample

	#Treatment
pred.t.ob=raw.ob%*%coeff.t
pred.t.in=raw.in%*%coeff.t

#This subtracts predicted values using matrix subtraction to get pites per bootstrap
pite.ob=pred.t.ob-pred.c.ob 
pite.in=pred.t.in-pred.c.in

#This is the mean pite across bootstraps
mpite.ob=apply(pite.ob,1,mean)
mpite.in=apply(pite.in,1,mean)

names(data)
dat.ob=data[which(data$ob==1),c(76:81)]
dat.in=as.data.frame(data[which(data$ob==0),c(76:81)])
#AGGREGATE PREDICTED EFFECTS#
dat.ob$mpite=mpite.ob
dat.in$mpite=mpite.in
dat.ob$pred.t=pred.t.ob
dat.ob$pred.c=pred.c.ob
dat.in$pred.t=pred.t.in
dat.in$pred.c=pred.c.in

names(dat.ob)
filename=rbind(dat.ob, dat.in)
#COLLECT ALL RESULTS#


 return(filename)
}
)

## Step 5: Remember to stop the cluster!
stopCluster(cl)

## Check that the parallel and serial output are the same
all.equal(serial.output, par.output)
## [1] TRUE


#time

str(par.output)
PredictedEF=do.call(rbind.data.frame,par.output) 
  #unlist all objects
PredictedEF$sim=rep(1:Nsims,each=dim(PredictedEF)[1]/Nsims) 
  #add simulation number on each.
names(PredictedEF)
pred.c=(PredictedEF$pred.c)
pred.t=(PredictedEF$pred.t)
Predicted=PredictedEF[,1:7]
fix(Predicted)
write.table(PredictedEF, "I:\\Differential TX effects\\PredictedEF.txt",row.names=F)
ob=PredictedEF[which( PredictedEF$ob==1),]


#preliminary bias runs
ob$bias=ob$TE-ob$mpite
mean(ob$bias)
a=by(ob,ob$sim,function(x){mean(x$bias)})
attributes(a)
hist(a)
summary(a)
names(ob)
hist(a, breaks=10, main="Correlation between true and predicted across sims", xlab="correlation")


###############ESTIMATION QUALITY
#RMSE

#rmse=var+ bias2
#var=rmse-bias2

str(ob)
#this subsets the data for ease
rmse.dat=ob[,c(2,6,7)]
rmse.dat$id=rep(1:5000,10)#adds id variable
rmse.dat$err=(rmse.dat$mpite - rmse.dat$TE)^2#deviation from true
names(rmse.dat)
rwide=rmse.dat[,c(2,4,5)] #this makes the apply easier (don't need to skipcolumns)
rmse.wide=reshape(rwide,timevar="sim", idvar=c("id"), direction="wide")
names(rmse.wide)
rmse.wide$mse=apply(rmse.wide[,2:101],1,(mean))
rmse.wide$rmse=sqrt(rmse.wide$mse)
hist(rmse.wide$rmse)
length(rmse.wide$rmse)

#plots in paper draft
rm <- ggplot(rmse.wide, aes(x = rmse))
rm +geom_density()+  labs(x = "RMSE", y="Density")+
theme_bw(base_size = 12)


#VARIABILITY
#Variability is the deviation from the mean, need mean mpite
names(ob)
names(rmse.wide)
var.dat=ob[,c(2,6,7)]
var.dat$id=rep(1:5000,10)#adds id variable
names(var.dat)
var.wide=reshape(var.dat,timevar="sim", idvar=c("id", "TE"), direction="wide")
var.wide$mean=apply(var.wide[,3:102],1,mean)
names(var.wide)
.2765898-.6043022
.5840406-.6043022
-.2635286--.5892625
-.5698415--.5892625
deviations2=apply(var.wide[,3:102],2,function(x){
dev2=(x-var.wide$mean)^2
v=as.data.frame(cbind(var.wide[1:2],dev2))
return(v)}
)
devs=as.data.frame(matrix(unlist(deviations2),nrow=5000,byrow=F))
dev.df=(devs[,c(1:2,seq(3,300,by=3))])
str(dev.df)
dev.df$variance=apply(dev.df[,3:102],1,mean)
#fix(dev.df)
hist(dev.df$variance)
mean(variance)

#PLOT IN PAPER DRAFT
v <- ggplot(dev.df, aes(x = variance))
v +geom_density()+  labs(x = "Variance of estimator", y="Density")+
theme_bw(base_size = 12)

###BIAS 
#This is the within subject bias.
#bias2
names(ob)
bias2.dat=ob[,c(2,6,10)]
bias2.dat$id=rep(1:5000,10)#adds id variable
bias2.wide=reshape(bias2.dat,timevar="sim", idvar=c("id", "TE"), direction="wide")
names(bias2.wide)
bias2.wide$meanbias=apply(bias2.wide[,3:102],1,(mean))
#fix(bias2.wide)
mean(bias2.wide$meanbias)
names(bias2.wide)
hist(bias2.wide$meanbias)
plot(bias2.wide$meanbias,bias2.wide$TE)
abline(lm(bias2.wide$TE~bias2.wide$meanbias), col="red")


#bias plot in paper draft
b <- ggplot(bias2.wide, aes(x = meanbias))
b
b +geom_density()+  labs(x = "Average bias within subjects", y="Density")+
theme_bw(base_size = 12)



#check
#var=rmse-bias2
names(dev.df)
var=dev.df[,c(1:2,103)]
fix(var)
ms=rmse.wide[,c(1,102:103)]
fix(ms)
b=bias2.wide[,c(1,2,103)]
qual=merge(var,ms,by.x="V2",by.y="id")
names(qual)
qual=merge(qual, b,by.x="V2",by.y="id")

qual$biassq=qual$meanbias^2
qual$var=qual$mse-qual$biassq
names(qual)
qual[1:10,c(3,9)]
#This checks!


###########################PREDICTION
sample(1:100, 1)#random number generator, chose 5
#INDIVIDUAL SIMS
sim1=Predicted[which(Predicted$ob==1 & Predicted$sim==1),]
sim2=Predicted[which(Predicted$ob==1 & Predicted$sim==2),]
sim5=Predicted[which(Predicted$ob==1 & Predicted$sim==5),]

names(sim1)
sim1[1:10,2:3]
sim2[1:10,2:3]
sim5[1:10,1:5]
dim(sim1)
sim2[1:10,1:5]
dim(sim1)


#TRUE VS PREDICTED (ONE SIM)
names(predicted.t)
predicted.t=data.frame(matrix(unlist(pred.t), nrow=1000000, byrow=F)) #this makes the output of EM a dataframe
predicted.c=data.frame(matrix(unlist(pred.c), nrow=1000000, byrow=F)) #this makes the output of EM a dataframe

predicted.t$mean.t=apply(predicted.t[1:100],1,mean)
predicted.t=cbind(predicted.t,Predicted)
predicted.c$mean.c=apply(predicted.c[1:100],1,mean)
predicted.c=cbind(predicted.c,Predicted)
sim1.pred.t=predicted.t[which(predicted.t$ob==1 & predicted.t$sim==1),]
sim1.pred.c=predicted.c[which(predicted.c$ob==1 & predicted.c$sim==1),]
sim5.pred.t=predicted.t[which(predicted.t$ob==1 & predicted.t$sim==5),]
sim5.pred.c=predicted.c[which(predicted.c$ob==1 & predicted.c$sim==5),]


#Scatterplots: True vs predicted under TREATMENT
spt=ggplot(data = sim5.pred.t,aes(y= sim5.pred.t$y.obs.tx, x = sim5.pred.t$mean.t ))
spt=spt+geom_point(size=2.5, color="Grey", shape=1)
spt=spt+ geom_smooth(method = "lm", se=FALSE, formula = y~x, size=1, size=1, color="black")+
theme_bw()  +
labs(x = "Predicted Y", y="Observed Y (treatment condition)")
#sp=sp+ geom_abline(intercept=0, slope=1, color="red", size=1.5)
spt

#Scatterplots True vs predicted under CONTROL
spc=ggplot(data = sim5.pred.c,aes(y= sim5.pred.c$y.obs.control, x = sim5.pred.c$mean.c ))
spc=spc+geom_point(size=2.5, color="Grey", shape=1)
spc=spc+ geom_smooth(method = "lm", se=FALSE, formula = y~x, size=1, size=1, color="black")+
theme_bw()  +
labs(x = "Predicted Y", y="Observed Y (treatment condition)")
#spc=spc+ geom_abline(intercept=0, slope=1, color="red", size=1.5)
spc
names(sim1.pred.c)
names(sim5)

#Scatterplots PITE
sppite=ggplot(data = sim5,aes(y= sim5$TE, x = sim5$mpite,colour = factor(sim5$obs.txt) ,shape = factor(sim5$obs.txt)))
sppite=sppite+geom_point()+labs(x="PITE (Predicted Individual Treatment Effect)", y="True Effect")
sppite=sppite+
scale_colour_manual(name = "Treatment Arm", breaks = c("0", "1"),
		labels = c("Control", "Treatment"),values = c(("blue"),("grey68"))) +
scale_shape_manual(name = "Treatment Arm",
				 breaks = c("0", "1"),
                      labels = c("Control", "Treatment"),
                     	values = c(0,1))
sppite=sp+geom_smooth(method = "lm", se=FALSE, formula = y~x, color="black", size=1) +theme_classic()
sppite



#plotting regression coefficients
f=by(ob,ob$sim,function(x) summary(lm(x$TE~x$mpite)))
(f)
f2=lapply(f,function(x) as.numeric(x$coefficients[2]))
f2=as.numeric(unlist(f2))
dataf2=data.frame(f2)
hist(f2)
reg <- ggplot(dataf2, aes(x = f2))
reg=reg +geom_density()+  labs(x = "Regression coefficients across sims", y="Density")+
theme_bw(base_size = 12)



#Examination of predicted values
names(sim1)
edit(pred.t)

#Scatterplots
sp=ggplot(data = sim1,aes(y= sim1$, x = sim1$y.obs.tx))
sp=sp+geom_point(size=2.5, color="Grey", shape=1)
sp=sp+ geom_smooth(method = "lm", se=FALSE, formula = y~x, size=1.5,color="Medium Blue", size=1)+theme_bw()        
sp=sp+ geom_abline(intercept=0, slope=1, color="red", size=1.5)
sp


###############################################################################
#FIGURE 1 (this code is copied from above, where separate plots were drawn
library(gridExtra)
library(grid)
b <- ggplot(bias2.wide, aes(x = meanbias))
b<-b +geom_density()+  labs(x = "Mean bias", y="Density")+
theme_classic(base_size = 10)

v <- ggplot(dev.df, aes(x = variance))
v<-v +geom_density()+  labs(x = "Variance of estimator", y="Density")+
theme_classic(base_size = 10)

rm <- ggplot(rmse.wide, aes(x = rmse))
rm<-rm +geom_density()+  labs(x = "RMSE", y="Density")+
theme_classic(base_size = 10)

grid.arrange(b, v, rm, ncol=2)






FIGURE 2

#Scatterplots: True vs predicted under TREATMENT
spt=ggplot(data = sim5.pred.t,aes(y= sim5.pred.t$y.obs.tx, x = sim5.pred.t$mean.t ))
spt=spt+geom_point(size=2.5, color="Grey", shape=1)
spt=spt+ geom_smooth(method = "lm", se=FALSE, formula = y~x, size=1, size=1, color="black")+
theme_bw()  +
labs(x = "Predicted Y", y="Observed Y (treatment condition)")
#sp=sp+ geom_abline(intercept=0, slope=1, color="red", size=1.5)
spt

#Scatterplots True vs predicted under CONTROL
spc=ggplot(data = sim5.pred.c,aes(y= sim5.pred.c$y.obs.control, x = sim5.pred.c$mean.c ))
spc=spc+geom_point(size=2.5, color="Grey", shape=1)
spc=spc+ geom_smooth(method = "lm", se=FALSE, formula = y~x, size=1, size=1, color="black")+
theme_bw()  +
labs(x = "Predicted Y", y="Observed Y (treatment condition)")
#spc=spc+ geom_abline(intercept=0, slope=1, color="red", size=1.5)
spc
names(sim1.pred.c)
names(sim5)

#Scatterplots PITE
sppite=ggplot(data = sim5,aes(y= sim5$TE, x = sim5$mpite,colour = factor(sim5$obs.txt) ,shape = factor(sim5$obs.txt)))
sppite=sppite+geom_point()+labs(x="PITE", y="True Effect")
sppite=sppite+
scale_colour_manual(name = "Treatment Arm", breaks = c("0", "1"),
		labels = c("Control", "Treatment"),values = c(("blue"),("grey68"))) +
scale_shape_manual(name = "Treatment Arm",
				 breaks = c("0", "1"),
                      labels = c("Control", "Treatment"),
                     	values = c(0,1))
sppite=sppite+geom_smooth(method = "lm", se=FALSE, formula = y~x, color="black", size=1) +theme_classic()
sppite

grid.arrange(spt, spc, sppite, ncol=2)

FIGURE 3 Regs
reg <- ggplot(dataf2, aes(x = f2))
reg=reg +geom_density()+  labs(x = "Regression coefficients across sims", y="Density")+
theme_bw(base_size = 12)
reg

