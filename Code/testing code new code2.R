#LOAD LIBRARIES#

library(ggplot2)
library(snow) #for parallelization
library(mi)

mean(y)=0

#The number of samples from the mixture distribution
N = 100000                 

#Sample N random uniforms U
U =runif(N)
hist(U)

#Variable to store the samples from the mixture distribution                                             
rand.samples = rep(NA,N)

#Sampling from the mixture
for(i in 1:N){
    if(U[i]<.9){
        rand.samples[i] = rnorm(1,0,.20)
        }else{
        rand.samples[i] = rnorm(1,-3,.5)
    }
}

logit

(1/(1+exp(-1*(15))))

#Density plot of the random samples
plot(density(rand.samples),main="Density Estimate of the Mixture Model")

#Plotting the true density as a sanity check
x = seq(-20,20,.1)
truth = .3*dnorm(x,0,1) + .5*dnorm(x,10,1) + .2*dnorm(x,3,.1)
plot(density(rand.samples),main="Density Estimate of the Mixture Model",ylim=c(0,.2),lwd=2)
lines(x,truth,col="red",lwd=2)

legend("topleft",c("True Density","Estimated Density"),col=c("red","black"),lwd=2)











dataGen=function(N,p,signp,seed){
 set.seed(2398)
Xsign=rbinom(N*signp,1,.05)
df1=data.frame(matrix(Xsign,nrow=N, ncol=signp))

Xns=rbinom(N*(p-signp),1,.5)
 df2=data.frame(matrix(Xns,nrow=N,ncol=(p-signp)))

df=cbind(df1,df2)
names(df) <- paste("X", 1:p, sep="")
 df$obs.txt=rep(0:1,N/2)
df$TE=(-4.5+1.1*df$X1+
	  1.4*df$X2+
	  2.5*df$X3+
	   3.3*df$X4+
         2.2*df$X5+
	  1.6*df$X6+
	  1.8*df$X7)
#These are logits(probabilities) of experiencing the event if you have the gene

#prob of experiencing event if no treatment (control) = .02
#prob of experiencing event if treatment = pTE
#How do I get treatment effects with negative direction???
	#df$oddsTE=exp(df$TE)
	df$pTE= (1/(1+exp(-1*(df$TE))))
	#df$pTEoriginal=df$oddsTE/(1+df$oddsTE)
 seed=set.seed(seed)
 df$y.obs.c=rbinom(N,1,.01)#observed y value under control
 df$y.obs.tx= ifelse(df$obs.txt==1, rbinom(N,1, df$pTE),NA) #observed y value under TX
 df$Y=ifelse(df$obs.txt==0,df$y.obs.c,df$y.obs.tx) #observed Y value
 df$y.obs.control=ifelse(df$obs.txt==0,df$y.obs.c,NA) #observed y value under control
 df$ob=rep(0:1,each=N/2)
 df$sim=rep(length(seed),each=N)
 return(df)}


#testing
a=dataGen(10000,76,7,543)
hist(a$pTE-.02)
hist(a$TE)
hist(a$Y)
mean(a$y.obs.tx, na.rm=T)
mean(a$y.obs.control, na.rm=T)
round(table(a$pTEneg),1)
x=a
names(x)
#CREATE MI FUNCTION#
mi.func<-function(x){
txt.imp=as.data.frame(x[which(x$ob==0),c(81,1:p)])
cont.imp=as.data.frame(x[which(x$ob==0),c(83,1:p)])
names(txt.imp)
names(cont.imp)

# STEP 1: CONVERT IT TO A missing_data.frame
mdf <- missing_data.frame(txt.imp)
mdf2 <- missing_data.frame(cont.imp)
	#show(mdf)
	#summary(mdf)
	#image(mdf)
	#str(mdf)

#This does the imputation for the control condition (imputes the treatment)
mi.control.i<-mi (mdf,  max.minutes=20000, n.chains=nimp)
mi.txt.i    <-mi (mdf2, max.minutes=20000, n.chains=nimp)

	#str(mi.control.i)

Conv.control=Rhats(mi.control.i)
Conv.txt=Rhats(mi.txt.i)


#gives complete data
comp.c=complete(mi.control.i,m=nimp)

tf <- vector("list", nimp)
for (j in 1:nimp){ 
        s=noquote(paste0("chain:",j))
		s
      tf[[j]] =  mi.control.i@data[[s]]@variables$y.obs.tx@parameters[30,]
}
tf
c=do.call(rbind,tf) #converts list to matrix

tf2 <- vector("list", nimp)
for (k in 1:nimp){ 
        s=noquote(paste0("chain:",k))
		s
      tf2[[k]] =  mi.txt.i@data[[s]]@variables$y.obs.control@parameters[30,]
}
tf2
b=do.call(rbind,tf2)
#Each row is the coefficients for that imputation.
parms=list(c,b,Conv.control, Conv.txt)
return(parms)
}


#RUN ANALYSES
Nsims=1
  #number of simulatins
N=10000
  #number of cases per sim
p=75
signp=7
#nsp=p-signp
  #number of predictors 
nimp=25
  #number of imputations per sim
predictors=c(1:75)
  #predictors that go into the model
set.seed(740)
seedlist=round(runif(Nsims,0,10000))

## Step 1: Create a cluster of child processes 
cl <- makeCluster(10)

## Step 2: Load the necessary R package(s)
## N.B. length(cl) is the number of child processes
##      in the cluster 

par.setup <- parLapply( cl, 1:length(cl),
    function(x) {
     	require(mi)})

## Step 3: Distribute the necessary R objects 
clusterExport( cl, c('seedlist','dataGen', 'mi.func','Nsims','N','p','nimp','predictors', 'signp'))

## Step 4: Do the computation
par.output <- parLapply(cl, seedlist,
    function(x) {
       filename=x

data=x
 data=dataGen(N,p,signp,x) #generate 100 cases with 75 predictors with seeds = x (aka "seedlist")
mi= mi.func(data)
str(mi)
coeff.t=(matrix(unlist(mi[1]),nrow=(p+1),byrow=T))
coeff.c=(matrix(unlist(mi[2]),nrow=(p+1),byrow=T))

#Raw data
raw.ob=as.matrix(data[which(data$ob==1),1:p])
raw.ob=as.matrix(cbind(rep(1,(dim(raw.ob)[1])),raw.ob))
raw.in=as.matrix(data[which(data$ob==0),1:p])
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

#This is the mean pite across imputations
mpite.ob=apply(pite.ob,1,mean)
mpite.in=apply(pite.in,1,mean)

#Converting logits to probabilites before subtracting
prob.pite.ob=(1/(1+exp(-1*(pred.t.ob))))-(1/(1+exp(-1*(pred.c.ob))))
prob.pite.in=(1/(1+exp(-1*(pred.t.in))))-(1/(1+exp(-1*(pred.c.in))))

mpite.prob.ob=apply(prob.pite.ob,1,mean)
mpite.prob.in=apply(prob.pite.in,1,mean)

mpite.prob.obneg=mpite.prob.ob*-1

names(data)
dat.ob=data[which(data$ob==1),c(((p+1)):85)]
dat.in=as.data.frame(data[which(data$ob==0),c((p+1):85)])

      #AGGREGATE PREDICTED EFFECTS#
dat.ob$mpite=mpite.ob
dat.in$mpite=mpite.in
dat.ob$pred.t=pred.t.ob
dat.ob$pred.c=pred.c.ob
dat.in$pred.t=pred.t.in
dat.in$pred.c=pred.c.in
dat.ob$mpite.prob=mpite.prob.ob
dat.in$mpite.prob=mpite.prob.in

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


plot(dat.ob$pTE,dat.ob$mpite.prob)
plot(dat.ob$TE,dat.ob$mpite)
hist(dat.ob$TE)
hist(dat.ob$mpite)
plot(dat.ob$pred.t.ob,dat.ob$mpite)
dim(pred.c.ob, dat.ob$y.ob.control)
hist(pred.c.ob)
#time

str(par.output)
PredictedEF=do.call(rbind.data.frame,par.output) 
  #unlist all objects
PredictedEF$sim=rep(1:Nsims,each=dim(PredictedEF)[1]/Nsims) 
  #add simulation number on each.
names(PredictedEF)
str(PredictedEF)
pred.c=(PredictedEF$pred.c)
pred.t=(PredictedEF$pred.t)
Predicted=PredictedEF[,c(1:10,13)]
#PredictedEF$logittoprob=1/(1+exp(-PredictedEF$mpite))
fix(Predicted)
write.table(PredictedEF, "I:\\Differential TX effects\\PredictedEF_categorical_R2small.txt",row.names=F)
ob=PredictedEF[which( PredictedEF$ob==1),]
dim(ob)

names(ob)
ob[1:10,c(2,10)]

#preliminary bias runs
ob$bias=ob$TE-ob$mpite
ob$prob.bias=ob$pTE-ob$mpite.prob
mean(ob$bias)
mean(ob$prob.bias)
a=by(ob,ob$sim,function(x){mean(x$bias)})
attributes(a)
hist(a)
summary(a)
names(ob)
#hist(a, breaks=10, main="Correlation between true and predicted across sims", xlab="correlation")


###############ESTIMATION QUALITY
#RMSE

#rmse=var+ bias2
#var=rmse-bias2

str(ob)
#this subsets the data for ease
rmse.dat=ob[,c(2,9,10)]
rmse.dat$id=rep(1:5000,Nsims)#adds id variable
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


#VARIABILITY###################################################################################
#Variability is the deviation from the mean, need mean mpite
names(ob)

#LOGIT SCALE
var.dat=ob[,c(2,9,10)]
var.dat$id=rep(1:5000,Nsims)#adds id variable
var.wide=reshape(var.dat,timevar="sim", idvar=c("id", "TE"), direction="wide")
names(var.wide)
var.wide$mean=apply(var.wide[,3:Nsims+2],1,mean)
deviations2=apply(var.wide[,3:102],2,function(x){
dev2=(x-var.wide$mean)^2
v=as.data.frame(cbind(var.wide[1:2],dev2))
return(v)}
)
devs=as.data.frame(matrix(unlist(deviations2),nrow=5000,byrow=F))
names(devs)
dev.df=(devs[,c(1:2,seq(3,300,by=3))])
dev.df$variance=apply(dev.df[,3:102],1,mean)
names(dev.df)
hist(dev.df$variance)
mean(dev.df$variance)
mean(1/(1+exp(-dev.df$variance)))


#PLOT IN PAPER DRAFT
v <- ggplot(dev.df, aes(x = variance))
v +geom_density()+  labs(x = "Variance of estimator", y="Density")+
theme_bw(base_size = 12)


#PROB SCALE
names(ob)
var.dat.p=ob[,c(4,9,13)]
var.dat.p$id=rep(1:5000,Nsims)#adds id variable
var.wide.p=reshape(var.dat.p,timevar="sim", idvar=c("id", "pTE"), direction="wide")
names(var.wide.p)
var.wide.p$mean=apply(var.wide.p[,3:Nsims+2],1,mean)
deviations2.p=apply(var.wide.p[,3:102],2,function(x){
dev2.p=(x-var.wide.p$mean)^2
v=as.data.frame(cbind(var.wide.p[1:2],dev2.p))
return(v)}
)
devs.p=as.data.frame(matrix(unlist(deviations2.p),nrow=5000,byrow=F))
names(devs.p)
dev.df.p=(devs.p[,c(1:2,seq(3,300,by=3))])
dev.df.p$variance=apply(dev.df.p[,3:102],1,mean)
names(dev.df.p)
hist(dev.df.p$variance)
mean(dev.df.p$variance)




###BIAS ###########################################################################
#This is the within subject bias.
#bias2
names(ob)
names(bias2.dat)
bias2.dat=ob[,c(2,10,15)]
bias2.dat$id=rep(1:5000,Nsims)#adds id variable
bias2.wide=reshape(bias2.dat,timevar="sim", idvar=c("id", "TE"), direction="wide")
names(bias2.wide)

bias2.wide$meanbias=apply(bias2.wide[,3:102],1,(mean))
mean(bias2.wide$meanbias)
hist(bias2.wide$meanbias)
plot(bias2.wide$meanbias,bias2.wide$TE)
abline(lm(bias2.wide$TE~bias2.wide$meanbias), col="red")

names(ob)
plot(ob$TE,ob$mpite)
abline(lm(ob$TE~ob$mpite), col="red")
plot(ob$pTE,ob$mpite.prob, xlim=c(-1,1),ylim=c(-1,1))
abline(lm(ob$pTE~ob$mpite.prob), col="red")


convertedbias=1/(1+exp(-1*(bias2.wide$meanbias)))
hist(convertedbias)
f=as.data.frame(cbind(convertedbias,ob$logittoprob))
names(f)
plot(f$convertedbias,f$V2)

###PROB SCALE
names(ob)
(#be sure to run prelim bias runs above)
bias2.dat.p=ob[,c(4,9,16)]
dim(ob)
names(bias2.dat.p)
bias2.dat.p$id=rep(1:5000,Nsims)#adds id variable
#bias2.dat.p$bias=bias2.dat.p$pTE-bias2.dat.p$mpite.prob
bias2.wide.p=reshape(bias2.dat.p,timevar="sim", idvar=c("id", "pTE"), direction="wide")
names(bias2.wide.p)

bias2.wide.p$meanbias.p=apply(bias2.wide.p[,3:102],1,(mean))

mean(bias2.wide.p$meanbias.p)
hist(bias2.wide.p$meanbias.p)
#fix(bias2.wide.p)
plot(bias2.wide.p$meanbias.p,bias2.wide.p$pTE)
abline(lm(bias2.wide.p$pTE~bias2.wide.p$meanbias.p), col="red")


#bias plot in paper draft
b <- ggplot(bias2.wide, aes(x = meanbias))
b +geom_density()+  labs(x = "Average bias within subjects", y="Density")+
theme_bw(base_size = 12)
names(bias2.wide.p)
b <- ggplot(bias2.wide.p, aes(x = meanbias.p))
b +geom_density()+  labs(x = "Average bias within subjects", y="Density")+
theme_bw(base_size = 12)

ob[1:10,c(4,13,16)]
ob[1:10,c(2,10,15)]

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
str(pred.t)
predicted.t$mean.t=apply(predicted.t[1:Nsims],1,mean)
predicted.t=cbind(predicted.t,Predicted)
predicted.c$mean.c=apply(predicted.c[1:Nsims],1,mean)
predicted.c=cbind(predicted.c,Predicted)
sim1.pred.t=predicted.t[which(predicted.t$ob==1 & predicted.t$sim==1),]
sim1.pred.c=predicted.c[which(predicted.c$ob==1 & predicted.c$sim==1),]
sim5.pred.t=predicted.t[which(predicted.t$ob==1 & predicted.t$sim==5),]
sim5.pred.c=predicted.c[which(predicted.c$ob==1 & predicted.c$sim==5),]


#Scatterplots: True vs predicted under TREATMENT
mean(sim5$y.obs.tx, na.rm=T)
mean(sim5.pred.t$mean.t)

mean(sim5$y.obs.c, na.rm=T)
mean(sim5.pred.c$mean.c)



#Scatterplots PITE
sppite=ggplot(data = sim1,aes(y= sim1$pTE, x = sim1$mpite.prob,colour = factor(sim1$obs.txt) ,shape = factor(sim1$obs.txt)))
sppite=sppite+geom_point()+labs(x="PITE (Predicted Individual Treatment Effect)", y="True Effect")
sppite=sppite+
scale_colour_manual(name = "Treatment Arm", breaks = c("0", "1"),
		labels = c("Control", "Treatment"),values = c(("blue"),("grey68"))) +
scale_shape_manual(name = "Treatment Arm",
				 breaks = c("0", "1"),
                      labels = c("Control", "Treatment"),
                     	values = c(0,1))
sppite=sppite+geom_smooth(method = "lm", se=FALSE, formula = y~x, color="black", size=1) +theme_classic()
sppite


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

#Scatterplots PITE
sppite=ggplot(data = sim1,aes(y= sim1$pTE-.02, x = sim1$mpite.prob,colour = factor(sim1$obs.txt) ,shape = factor(sim1$obs.txt)))
sppite=sppite+geom_point()+labs(x="PITE (Predicted Individual Treatment Effect),IMP", y="True Effect")
sppite=sppite+
scale_colour_manual(name = "Treatment Arm", breaks = c("0", "1"),
		labels = c("Control", "Treatment"),values = c(("dodgerblue2"),("yellow3"))) +
scale_shape_manual(name = "Treatment Arm",
				 breaks = c("0", "1"),
                      labels = c("Control", "Treatment"),
                     	values = c(0,1))
sppite=sppite+geom_smooth(method = "lm", se=FALSE, formula = y~x,  size=1) +theme_classic()+ylim(-1,1)+xlim(-1,1)+
geom_abline(slope=1, intercept=0, col="red")
sppite

grid.arrange(spt, spc, sppite, ncol=2)
