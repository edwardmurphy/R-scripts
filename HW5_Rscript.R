setwd("C://Users//Owner//Documents//STAT 676")
Sweave("HW5_Sweave.nw")

library(xtable)
library(arm)
library(BRugs)

x<-seq(-10,10,by=0.001)
mixden<-function(x){
  exp(-2*(x+1)^2) + exp(-0.5*(x-1)^2)
}
plot(x,mixden,type="l",ylim=c(0,1.5))
lines(x,dnorm(x),lty=2)
lines(x,dnorm(x,sd=2),lty=3)
lines(x,7*dnorm(x,sd=2),lty=4)

xstar<-runif(10000,-5,5)
sdstar<-runif(10000,.5,2)
maxM<-rep(0,times=length(sdstar))
for(i in 1:length(sdstar)){
  m<-mixden(xstar)/dnorm(xstar,sd=sdstar[[i]])
  maxM[[i]]<-max(m)
}

sdandm<-cbind(sdstar,maxM)
sdandm[which.min(sdandm[,2]),]

foverg<-function(x,sig){
  (exp(-2*(x+1)^2) + exp(-0.5*(x-1)^2))/dnorm(x,sd=sig)
}


##I think this is wrong b/c we want to find the maximin M, so maybe separate into two problems (one for sd, other for M)
c<-5; rSD<-0.5; rX<-1
T<-c; sd<-1; x<-0
iter<-10000
for(j in 1:iter){
  a<-max(sd-rSD,0.001)
  uSD<-runif(1,min=a,max=sd+rSD)
  uX<-runif(1,min=x-rX,max=x+rX)
  rho<-exp(-(foverg(uX,uSD)-foverg(x,sd))/T)
  v<-runif(1)
  sd<-if(v<=rho)uSD else sd
  x<-if(v<=rho)uX else x
  T<-c/log(j+1)
}







