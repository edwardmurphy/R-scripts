setwd("C://Users//Owner//Documents//STAT 676")
Sweave("HW3_Sweave.nw")

library(xtable)

###Problem 1

a<-seq(2,20,by=1)
b<-14.5/(a-1)

paramSearch<-cbind(a,b,qgamma(p=0.90,a,scale=b))
colnames(paramSearch)<-c("alpha","beta","90th Percentile")

###Problem 2

gamdes<-function(alpha,beta){
mean<-alpha*beta
mode<-beta*(alpha-1)
credSet<-qgamma(rbind(0.025,0.975),alpha,scale=beta)
output<-rbind(mean,mode,credSet)
rownames(output)<-c("mean","mode","2.5th Percentile","97.5th Percentile")
return(output)
}

prior<-gamdes(15,1)
post<-gamdes(97,1/8)

distSummary<-cbind(prior[1:2,],post[1:2,])
colnames(distSummary)<-c("Prior","Posterior")

xtable(distSummary,digits=2)

###Problem 3

post<-gamdes(82.5,1/7)
colnames(post)<-c("Posterior")

xtable(post)

###Problem 4

post<-gamdes(83,1/7)
colnames(post)<-c("Posterior")

xtable(post)

###Problem 5

lam<-seq(0.001,30,by=0.001)

prior1<-dgamma(lam,15,scale=1)
post1<-dgamma(lam,97,scale=1/8)

prior2<-sqrt(1/lam)
post2<-dgamma(lam,82.5,scale=1/7)

prior3<-rep(0.5,length(lam))
post3<-dgamma(lam,83,scale=1/7)

par(mfrow=c(1,3))
plot(lam,post1,type="l",ylim=c(0,0.5),xlab="Lambda",ylab="Density",main="Gamma(15,1) Prior")
lines(lam,prior1,lty=2)

plot(lam,post2,type="l",ylim=c(0,0.5),xlab="Lambda",ylab="Density",main="Jeffreys Prior")
lines(lam,prior2,lty=2)

plot(lam,post3,type="l",ylim=c(0,0.5),xlab="Lambda",ylab="Density",main="Uniform Prior")
lines(lam,prior3,lty=2)

par(mfrow=c(1,1))
plot(lam,post1,type="l",xlab="Lambda",ylab="Density")
lines(lam,post2,lty=2,col="blue")
lines(lam,post3,lty=5,col="red")
legend(16,0.3,c("Gamma(15,1) Prior","Jeffreys Prior","Uniform Prior"),
lty=c(1,2,5),col=c("black","blue","red"))