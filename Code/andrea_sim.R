#Andrea SIM#
sim.andrea.data=function(N=10000,p=7){
  
  
  #binary covariates generated to resemble design of the data in the motivating example
  set.seed(2398)
  X=rbinom(N*p,1,.5)
  
  df.c=matrix(X,nrow=N,ncol=p)
  df.d=matrix(0,nrow=N,ncol=68)
  df=data.frame(df.c,df.d)
  colnames(df) <- paste("X",1:ncol(df),sep="")
  
  # split half-half, 50% for derivation of predictive model, 50% for out-of-sample calculation
  df$obs.txt=rep(0:1,N/2)
  
  # True treatment effect of each individual 
  df$TE=-1.3*df$X1-1.2*df$X2-.6*df$X3+.3*df$X4+.5*df$X5+1.1*df$X6+1.2*df$X7
  
  #observed y value under control
  set.seed(2984)
  df$y.obs.control=rnorm(N,0,1)
  
  #observed y value under TX
  df$y.obs.tx= ifelse(df$obs.txt==1, (df$y.obs.control + df$TE),NA) 
  
  #observed Y value
  df$Y=ifelse(df$obs.txt==0,df$y.obs.control,df$y.obs.tx) 
  
  #observed y value under control
  df$y.obs.control=ifelse(df$obs.txt==0,df$y.obs.control,NA) 
  
  df$ob=rep(0:1,each=10/2)
  return(df)
}