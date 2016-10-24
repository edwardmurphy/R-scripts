###################################################
### chunk number 1: 
###################################################
#line 28 "HW3_Sweave.nw"
library(xtable)


###################################################
### chunk number 2: 
###################################################
#line 47 "HW3_Sweave.nw"
a<-seq(2,20,by=1)
b<-14.5/(a-1)

paramSearch<-cbind(a,b,qgamma(p=0.90,a,scale=b))
colnames(paramSearch)<-c("alpha","beta","90th Percentile")


###################################################
### chunk number 3: 
###################################################
#line 55 "HW3_Sweave.nw"
xtable(paramSearch[9:19,],digits=1,caption='Grid Search Results for Gamma Hyperparameters',label="tab:1a")


###################################################
### chunk number 4: 
###################################################
#line 76 "HW3_Sweave.nw"
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


###################################################
### chunk number 5: 
###################################################
#line 93 "HW3_Sweave.nw"
xtable(distSummary,digits=2,caption='Prior and Posterior Statistics',label="tab:1c")


###################################################
### chunk number 6: 
###################################################
#line 127 "HW3_Sweave.nw"
post<-gamdes(82.5,1/7)
colnames(post)<-c("Posterior")


###################################################
### chunk number 7: 
###################################################
#line 132 "HW3_Sweave.nw"
xtable(post,digits=2,caption="Posterior Mean, Mode and 95 Percent Credible Set from Jeffreys Prior",label="tab:2c")


###################################################
### chunk number 8: 
###################################################
#line 150 "HW3_Sweave.nw"
post<-gamdes(83,1/7)
colnames(post)<-c("Posterior")


###################################################
### chunk number 9: 
###################################################
#line 155 "HW3_Sweave.nw"
xtable(post,digits=2,caption='Posterior Mean, Mode and 95 Percent Credible Set from Uniform Prior',label="tab:3b")


###################################################
### chunk number 10: fig5a
###################################################
#line 163 "HW3_Sweave.nw"
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


###################################################
### chunk number 11: 5a
###################################################
#line 187 "HW3_Sweave.nw"
#line 163 "HW3_Sweave.nw#from line#187#"
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
#line 188 "HW3_Sweave.nw"


###################################################
### chunk number 12: fig5b
###################################################
#line 195 "HW3_Sweave.nw"
par(mfrow=c(1,1))
plot(lam,post1,type="l")
lines(lam,post2,lty=2,col="blue")
lines(lam,post3,lty=5,col="red")
legend(16,0.3,c("Gamma(15,1) Prior","Jeffreys Prior","Uniform Prior"),
lty=c(1,2,5),col=c("black","blue","red"))


###################################################
### chunk number 13: 5b
###################################################
#line 205 "HW3_Sweave.nw"
#line 195 "HW3_Sweave.nw#from line#205#"
par(mfrow=c(1,1))
plot(lam,post1,type="l")
lines(lam,post2,lty=2,col="blue")
lines(lam,post3,lty=5,col="red")
legend(16,0.3,c("Gamma(15,1) Prior","Jeffreys Prior","Uniform Prior"),
lty=c(1,2,5),col=c("black","blue","red"))
#line 206 "HW3_Sweave.nw"


