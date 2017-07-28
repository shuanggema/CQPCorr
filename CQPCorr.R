##########################################################################
## Censored Quantile Partial Correlation (CQPCorr) 

## Publication: 
## Robust Identification of Gene-environment Interactions for 
## Prognosis Using a Censored Quantile Partial Correlation Approach
## by Yaqing Xu, Mengyun Wu, Qingzhao Zhang, and Shuangge Ma (2017)

## Maintainer: Yaqing Xu <yaqing.xu@yale.edu>
##########################################################################

############ 1. Functions ############
reweight<-function(y,delta,tau){
  # calculate weights for step I in cqpcorr
  # y is the survival time
  # delta is the censoring status
  # tau is the quantile level of interest
  kmtau<-survfit(Surv(y,delta)~1,type="kaplan-meier")
  fhat<-(1-kmtau$surv)[order(order(y))]

  w<-rep(1,n)
  index<-which(delta==0) 
  if(length(index)>0){
    for(i in 1:length(index)){
      if(fhat[index[i]]<tau) w[index[i]]<-(tau-fhat[index[i]])/(1-fhat[index[i]])
    }}
  return(w)
}


cqpcorr<-function(y,delta,x,z,w,tau){
  # note that x and z are univariates 
  # w is the weights for step I
  # x and z are representing environmental and genetic factors respectively
  
  inter<-x*z

  ### step I
  #construct pseudo observations for y^\infty
  index<-which(w!=1)
  y.pse<-rep(max(y)+9999,length(index))
  x.pse<-x[index]
  z.pse<-z[index]
  
  yy<-c(y,y.pse)
  xx<-c(x,x.pse)
  zz<-c(z,z.pse)
  ww<-c(w,1-w[index])
  
  part2<-rq(yy~xx+zz+1,weights=ww,tau)
  
  ### step II
  part1<-lm(inter~x+z+1)
  
  ### step III: quantile partial corr
  ind<-as.numeric(part2$residuals[1:n]<0)
  var1<-tau-w*ind
  var2<-part1$residuals
  yxz<-(mean(var1*var2)-mean(var1)*mean(var2))/sqrt((mean(w^2)*tau-mean(w)^2*tau^2)*var(var2))
  
  return(yxz)
}


################ 2.Example ################
library(survival)
library(quantreg)
library(MASS)
p=5 
q=3
n=200

z<-mvrnorm(n,rep(0,p),diag(p))
x<-mvrnorm(n,rep(0,q),diag(q))
y<-x%*%rep(1,q)+z%*%rep(2,p)+rnorm(n)
delta<-rep(1,n)
tau<-0.5

cqpcorr_results<-matrix(0,nrow=p,ncol=q)

w<-reweight(y,delta,tau)
for(m in 1:q){
  for(k in 1:p){
    cqpcorr_results[k,m]<-cqpcorr(y,delta,x[,m],z[,k],w,tau)
  }
}
cqpcorr_results
