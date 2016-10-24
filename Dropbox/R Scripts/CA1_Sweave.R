###################################################
### chunk number 1: 
###################################################
library(xtable)
x<-c(0:12)


###################################################
### chunk number 2: fig1ab
###################################################
par(mfrow=c(1,2))
plot(x,dbinom(x,12,0.25),type="h",lwd=2,ylab="density",main="pmf")
plot(x,pbinom(x,12,0.25),type="l",lwd=2,ylab="cumulative density",main="cdf")


###################################################
### chunk number 3: fig1ab
###################################################
par(mfrow=c(1,2))
plot(x,dbinom(x,12,0.25),type="h",lwd=2,ylab="density",main="pmf")
plot(x,pbinom(x,12,0.25),type="l",lwd=2,ylab="cumulative density",main="cdf")


###################################################
### chunk number 4: 
###################################################
rb<-rbinom(50,12,0.25)
rbHist<-hist(rb,breaks=seq(-0.5,12.5,by=1))


###################################################
### chunk number 5: fig1c
###################################################
par(mfrow=c(1,1))
plot(rbHist,xlab="x",main="",freq=F)


###################################################
### chunk number 6: fig1c
###################################################
par(mfrow=c(1,1))
plot(rbHist,xlab="x",main="",freq=F)


###################################################
### chunk number 7: 
###################################################
empTab<-cbind(rbHist$density,dbinom(x,12,0.25))
colnames(empTab)<-c("empirical","true")
rownames(empTab)<-seq(0,12,by=1)
options(digits=4)


###################################################
### chunk number 8: 
###################################################
xtable(empTab,digits=4,caption="Empirical density of 50 random variates and true density for Binomial(12,0.25) distribution")


###################################################
### chunk number 9: 
###################################################
x<-seq(0,1,by=0.001)


###################################################
### chunk number 10: fig2a
###################################################
par(mfrow=c(1,2))
plot(x,dbeta(x,2,1),type="l",ylab="density",main="pdf")
plot(x,pbeta(x,2,1),type="l",ylab="cumulative density",main="cdf")


###################################################
### chunk number 11: fig2a
###################################################
par(mfrow=c(1,2))
plot(x,dbeta(x,2,1),type="l",ylab="density",main="pdf")
plot(x,pbeta(x,2,1),type="l",ylab="cumulative density",main="cdf")


###################################################
### chunk number 12: 
###################################################
pcrtBeta<-cbind(round(qbeta(0.025,2,1),4),round(qbeta(0.975,2,1),4))
colnames(pcrtBeta)<-c("2.5","97.5")
rownames(pcrtBeta)<-"percentile value"


###################################################
### chunk number 13: 
###################################################
pcrtBeta


###################################################
### chunk number 14: 
###################################################
rBeta<-rbeta(100,2,1)


###################################################
### chunk number 15: fig2d
###################################################
par(mfrow=c(1,1))
plot(density(rBeta),xlim=c(0,1),ylim=c(0,2),xlab="x",main="")
lines(x,dbeta(x,2,1),lty=2)
legend(0,2,c("kernel density","true density"),lty=1:2)


###################################################
### chunk number 16: fig2d
###################################################
par(mfrow=c(1,1))
plot(density(rBeta),xlim=c(0,1),ylim=c(0,2),xlab="x",main="")
lines(x,dbeta(x,2,1),lty=2)
legend(0,2,c("kernel density","true density"),lty=1:2)


###################################################
### chunk number 17: 
###################################################
y<-6
n<-8
darciPost<-dbeta(x,1+y,1+n-y)
ethanPost<-dbeta(x,4+y,4+n-y)
frankPost<-dbeta(x,6+y,2+n-y)


###################################################
### chunk number 18: fig3a
###################################################
par(mfrow=c(1,1))
plot(x,darciPost,type="l",ylim=c(0,max(frankPost)),ylab="density",main="Posterior Distribution")
lines(x,ethanPost,lty=2)
lines(x,frankPost,lty=3)

legend(0,3.5,c("Darci","Ethan","Frank"),lty=1:3)


###################################################
### chunk number 19: fig3a
###################################################
par(mfrow=c(1,1))
plot(x,darciPost,type="l",ylim=c(0,max(frankPost)),ylab="density",main="Posterior Distribution")
lines(x,ethanPost,lty=2)
lines(x,frankPost,lty=3)

legend(0,3.5,c("Darci","Ethan","Frank"),lty=1:3)


###################################################
### chunk number 20: 
###################################################
y<-15
n<-20

darciPrior<-dbeta(x,1,1)
ethanPrior<-dbeta(x,4,4)
frankPrior<-dbeta(x,6,2)

darciPost<-dbeta(x,1+y,1+n-y)
ethanPost<-dbeta(x,4+y,4+n-y)
frankPost<-dbeta(x,6+y,2+n-y)


###################################################
### chunk number 21: fig3b
###################################################
par(mfrow=c(1,3))
plot(x,darciPost,type="l",lty=2,ylab="density",main="Darci",ylim=c(0,5))
lines(x,darciPrior)

plot(x,ethanPost,type="l",lty=2,ylab="density",main="Ethan",ylim=c(0,5))
lines(x,ethanPrior)

plot(x,frankPost,type="l",lty=2,ylab="density",main="Frank",ylim=c(0,5))
lines(x,frankPrior)


###################################################
### chunk number 22: fig3b
###################################################
par(mfrow=c(1,3))
plot(x,darciPost,type="l",lty=2,ylab="density",main="Darci",ylim=c(0,5))
lines(x,darciPrior)

plot(x,ethanPost,type="l",lty=2,ylab="density",main="Ethan",ylim=c(0,5))
lines(x,ethanPrior)

plot(x,frankPost,type="l",lty=2,ylab="density",main="Frank",ylim=c(0,5))
lines(x,frankPrior)


###################################################
### chunk number 23: 
###################################################
betades<-function(shape1, shape2)
{
meanbeta<-shape1/(shape1+shape2)
modebeta<-(shape1-1)/(shape1+shape2-2)
medbeta<-qbeta(.5,shape1,shape2)
sdbeta<-sqrt(shape1*shape2/((shape1+shape2)^2*(shape1+shape2+1)))
des<-c(meanbeta, modebeta, medbeta, sdbeta)
names(des)<-c("mean", "mode", "median", "standard deviation")
return(des)
}

Darci<-betades(7,3)
Ethan<-betades(10,6)
Frank<-betades(12,4)


###################################################
### chunk number 24: 
###################################################
xtable(cbind(Darci,Ethan,Frank),digits=3,caption="Descriptive Statistics for Posterior Distributions from Problem 3a")


###################################################
### chunk number 25: 
###################################################
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

darciDes2<-betades2(7,3)
ethanDes2<-betades2(10,6)
frankDes2<-betades2(12,4)

des2.a<-cbind(darciDes2,ethanDes2,frankDes2)
colnames(des2.a)<-c("Darci","Ethan","Frank")


###################################################
### chunk number 26: 
###################################################
xtable(des2.a[5:6,],caption="2.5th and 97.5th Percentiles for Posterior Distributions from Problem 3a")


###################################################
### chunk number 27: 
###################################################

darciDes2.b<-betades2(16,6)
ethanDes2.b<-betades2(19,9)
frankDes2.b<-betades2(21,7)

des2.b<-cbind(darciDes2.b,ethanDes2.b,frankDes2.b)
colnames(des2.b)<-c("Darci","Ethan","Frank")


###################################################
### chunk number 28: 
###################################################
xtable(des2.b,caption="Descriptive Statistics for Posterior Ditributions in Problem 3b")


