setwd("C://Users//Owner//Documents//STAT676")
Sweave("HW2_Sweave.nw")

####Problem 1
library(xtable)
x<-seq(0,1,by=0.001)

a1 = 78; b1 = 12
a2 = 31; b2 = 3  ##gives prior mode of 94%
	           ##prior sample size is only 32, but 88 in prior1
a3 = 1; b3 = 1

prior1<-dbeta(x,a1,b1)
prior2<-dbeta(x,a2,b2)
prior3<-dbeta(x,a3,b3)

y<-22
n<-28

post1<-dbeta(x,a1+y,b1+n-y)
post2<-dbeta(x,a2+y,b2+n-y)
post3<-dbeta(x,a3+y,b3+n-y)

par(mfrow=c(1,3))
plot(x,post1,type="l",ylab="density",main="Part a",ylim=c(0,12))
lines(x,prior1,lty=2)

plot(x,post2,type="l",ylab="density",main="Part b",ylim=c(0,12))
lines(x,prior1,lty=2)

plot(x,post3,type="l",ylab="density",main="Part c",ylim=c(0,12))
lines(x,prior3,lty=2)

betades2<-function(shape1, shape2)
{
meanbeta<-shape1/(shape1+shape2)
modebeta<-(shape1-1)/(shape1+shape2-2)
medbeta<-qbeta(.5,shape1,shape2)
sdbeta<-sqrt(shape1*shape2/((shape1+shape2)^2*(shape1+shape2+1)))
prctile_2.5<-qbeta(0.025,shape1,shape2)
prctile_97.5<-qbeta(0.975,shape1,shape2)
des<-c(meanbeta, modebeta, medbeta, sdbeta, prctile_2.5, prctile_97.5)
names(des)<-c("mean", "mode", "median", "standard deviation","2.5th percentile","97.5th percentile")
return(des)
}

post1Des<-betades2(a1+y,b1+n-y)
post2Des<-betades2(a2+y,b2+n-y)
post3Des<-betades2(a3+y,b3+n-y)

postDes<-cbind(post1Des,post2Des,post3Des)
colnames(postDes)<-c("Part a","Part b","Part c")
xtable(postDes[c(1:2,5:6),],digits=3,caption="Mean and Mode for Posterior Distributions"))

postProb1<-1-pbeta(0.8,a1+y,b1+n-y)
postProb2<-1-pbeta(0.8,a2+y,b2+n-y)
postProb3<-1-pbeta(0.8,a3+y,b3+n-y)

postProb<-cbind(postProb1,postProb2,postProb3)
colnames(postProb)<-c("Part a","Part b","Part c")
rownames(postProb)<-c("P > 0.8")
xtable(postProb,digits=3,caption='Posterior Probability $\theta$ > 0.8')

####Problem 2
a = 6; b = 1

y<-69
n<-75

par(mfrow=c(1,1))
prior<-dbeta(x,a,b)
post<-dbeta(x,a+y,b+n-y)

plot(x,post,type="l",ylab="density",main="",ylim=c(0,14))
lines(x,prior,lty=2)

postDes<-cbind(betades2(a+y,b+n-y))
colnames(postDes)<-c("")
xtable(postDes,digits=3,caption="Posterior Statistics",label="tab:2d")

postOdds<-(1-pbeta(0.85,a+y,b+n-y))/pbeta(0.85,a+y,b+n-y)
priOdds<-(1-pbeta(0.85,6,1))/pbeta(0.85,6,1)
BF<-postOdds/priOdds

####Problem 3 (normal/normal, variance known)
priorMean<-0
priorVar<-2
dataVar<-2

postMean<-(priorMean/priorVar + y/dataVar)/(1/priorVar + 1/dataVar)
postVar<-1/(1/priorVar + 1/dataVar)

x<-seq(-10,10,by=0.001)
y<-4

prior<-dnorm(x,priorMean,sqrt(priorVar))
like<-dnorm(x,y,sqrt(dataVar))
post<-dnorm(x,postMean,sqrt(postVar))

par(mfrow=c(1,1))
plot(x,post,type="l",ylab="density")
lines(x,prior,lty=2,col="blue")
lines(x,like,lty=5,col="red")
legend(-10,0.4,c("Posterior","Prior","Likelihood"),lty=c(1,2,5),col=c("black","blue","red"))

##part b
priorVar<-18
postMean<-(priorMean/priorVar + y/dataVar)/(1/priorVar + 1/dataVar)
postVar<-1/(1/priorVar + 1/dataVar)

prior<-dnorm(x,priorMean,sqrt(priorVar))
like<-dnorm(x,y,sqrt(dataVar))
post<-dnorm(x,postMean,sqrt(postVar))

x11()
plot(x,post,type="l",ylab="density")
lines(x,prior,lty=2,col="blue")
lines(x,like,lty=5,col="red")
legend(-10,0.3,c("Posterior","Prior","Likelihood"),lty=c(1,2,5),col=c("black","blue","red"))


####Problem 4

#first let's generate some data using the prior hierarchy

nu0<-30
sig0<-4

#generate 30 sig2 
sig2<-1/rgamma(30,nu0/2,sig0*nu0/2)

#generate theta
theta<- rnorm(length(sig2),99,sqrt(sig2))

#generate data
y<-rnorm(length(theta),theta,sqrt(sig2))

hist(sig2)
x11()
hist(theta)
x11()
hist(y)

##now sample from posteriors
y<-rnorm(30,99,2) ##e.g. recovery data observed

nu0<-1
kap0<-1
mu0<-100 ##assume method is unbiased
sig0<-1  ##assume method is precise
n<-length(y)

nuN<-nu0+n
kapN<-kap0+n

sigN<-1/nuN*(nu0*sig0+(n-1)*var(y)+(kap0*n/kapN)*(mean(y)-mu0)^2)
sigN

#sig2 post
sig2<-1/rgamma(10000,nuN/2,nuN*sigN/2)

#theta post
muN<-(kap0*mu0 + n*mean(y))/(kapN)
theta<-rnorm(10000,muN,sqrt(sig2/kapN))
quantile(theta,0.025)
quantile(theta,0.975)

#t intervals
muN-qt(0.975,nuN)*sqrt(sigN/kapN)
muN+qt(0.975,nuN)*sqrt(sigN/kapN)
###very close!

#plot kernel
plot(density(theta),type="l")










