indx=NULL

for(i in 1:length(test.cntrl.data[,1])){
  indx=c(indx,sum(test.cntrl.data[i,1:7])==7)
}

quantile(test.cntrl.data[indx,8],na.rm=T)


vif_data <- estimator.cntrl.data[,1:7]
View(vif_data)
vif(vif_data)
vifstep(vif_data)
cor(vif_data)
mean((test.cntrl.data$control-test.cntrl.data$control2)^2,na.rm=T)/sd(test.cntrl.data$control,na.rm=TRUE)
mean((test.cntrl.data$control-test.cntrl.data$control2)^2,na.rm=T)/sd(test.cntrl.data$control,na.rm=TRUE)
summary(lm(control~.,data=test.cntrl.data))
summary(rfsrc(control~.,data=test.cntrl.data))
rfsrc(control~.,data=test.cntrl.data)
princomp
princomp(test.cntrl.data[,1:7])
pca=princomp(test.cntrl.data[,1:7])
summary(pca)
sum(test.cntrl.data[,1]%*%test.cntrl.data[,2])
sum(test.cntrl.data[,1]*test.cntrl.data[,2])
mean(test.cntrl.data[,1]*test.cntrl.data[,2])
mean(test.cntrl.data[,1]*test.cntrl.data[,5])
10000/128
which(test.cntrl.data[,1:7]==test.cntrl.data[,1:7])
which(test.cntrl.data[,1:7]==test.cntrl.data[1,1:7])
which(test.cntrl.data[,1:7]==1)
which(test.cntrl.data[,1:7]==rep(1,7))
test.cntrl.data[test.cntrl.data[,1:7]==rep(1,7),]

test.cntrl.data[sum(test.cntrl.data[,1:7])==7,]