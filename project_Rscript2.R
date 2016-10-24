### We will compare 2 models here,
## All models have uniform distribution on the mean slope and intercept parameters
## Model 1 has gamma dists on the precision (tau) parameters.  As Gelman has shown
## the use of this prior can heavily influence the posterior depending on the hyper
## parameters.  we will use (0.001,0.001) 

## Model 2 has uniform dists on the standard deviation parameters. We use what we believe are
## realistic values of the standard deviation based on the range 

setwd("C://Users//Owner//Documents//STAT 676")

library(xtable)
library(arm)
library(BRugs)
library(boa)

Sweave("presentation_Sweave.nw")

##Model 1a gamma(0.001,0.001)

inits1<-list(tau.g = 1 , tau.beta0 = 0.1, tau.beta1 = 0.1)#, mu.beta0 = 101, mu.beta1 = -0.3)
inits2<-list(tau.g = 1 , tau.beta0 = 1, tau.beta1 = 1)#, mu.beta0 = 101, mu.beta1 = -0.3)
inits3<-list(tau.g = 1 , tau.beta0 = 100, tau.beta1 = 100)#, mu.beta0 = 101, mu.beta1 = -0.3)
inits<-list(inits1,inits2,inits3)

parameters<-c("beta0[]","beta1[]","mu.beta0","mu.beta1","tau.g","tau.beta0","tau.beta1",
"sig2.g","sig2.beta0","sig2.beta1")

sim1a<-bugs(data="project_data.txt",inits,parameters,"project_bugsmodel.txt",
  debug=T,
  n.chains=3,n.iter=16000,n.burnin=1000,n.thin=1,
  bugs.seed=floor(runif(1,1,10000)))

sims1a<-sim1a$sims.matrix

##Model 1b gamma(1,1)

# sim1b<-bugs(data="project_data.txt",inits,parameters,"project_bugsmodel1b.txt",
#   #debug=T,
#   n.chains=3,n.iter=16000,n.burnin=1000,n.thin=1,
#   bugs.seed=floor(runif(1,1,10000)))
# 
# sims1b<-sim1b$sims.matrix

##Model 4a narrow uniforms
inits1<-list(sig.g = 0.1 , sig.beta0 = 0.1, sig.beta1 = 0.01)#, mu.beta0 = 101, mu.beta1 = -0.3)
inits2<-list(sig.g = 1 , sig.beta0 = 1, sig.beta1 = 0.1)#, mu.beta0 = 101, mu.beta1 = -0.3)
inits3<-list(sig.g = 10 , sig.beta0 = 5, sig.beta1 = 1)#, mu.beta0 = 101, mu.beta1 = -0.3)
inits<-list(inits1,inits2,inits3)

sim4a<-bugs(data="project_data.txt",inits,parameters,"project_bugsmodel4.txt",
  debug=T,
  n.chains=3,n.iter=16000,n.burnin=1000,n.thin=1,
  bugs.seed=floor(runif(1,1,10000)))

sims4a<-sim4a$sims.matrix

# ##Model 4b wider uniforms
# inits1<-list(sig.g = 0.1 , sig.beta0 = 0.1, sig.beta1 = 0.1)#, mu.beta0 = 101, mu.beta1 = -0.3)
# inits2<-list(sig.g = 1 , sig.beta0 = 1, sig.beta1 = 1)#, mu.beta0 = 101, mu.beta1 = -0.3)
# inits3<-list(sig.g = 100 , sig.beta0 = 50, sig.beta1 = 10)#, mu.beta0 = 101, mu.beta1 = -0.3)
# inits<-list(inits1,inits2,inits3)
# 
# sim4b<-bugs(data="project_data.txt",inits,parameters,"project_bugsmodel4b.txt",
#   #debug=T,
#   n.chains=3,n.iter=16000,n.burnin=1000,n.thin=1,
#   bugs.seed=floor(runif(1,1,10000)))
# 
# sims4b<-sim4b$sims.matrix

xtable(sim1a$summary[,c(1:3,5,7)],digits=4)

xtable(sim4a$summary[,c(1:3,5,7)],digits=4)


x11()
plot(density(sims1a[,18]),xlim=c(0,1),main="Posterior Density of Batch Slope Variance",
  xlab="")
lines(density(sims1b[,18]),lty=2)
lines(density(sims4a[,18]),lty=3)
lines(density(sims4b[,18]),lty=4)
legend("topright",legend=c("Model 1a","Model 1b","Model 4a",
  "Model 4b"),lty=1:4)

##batch slope density plots
x11()
par(mfrow=c(1,1))
#Model 1
plot(density(sims4a[,12]),lwd=3,xlim=c(-1,1),
  #main="Pooled and Individual Density Plots for Slope, Model 2",
  main="",
  xlab="Posterior Slope",ylim=c(0,7))
lines(density(sims4a[,6]),col="purple")
lines(density(sims4a[,7]),col="red")
lines(density(sims4a[,8]),col="blue",lty=2)
lines(density(sims4a[,9]),col="green")
lines(density(sims4a[,10]),col="orange",lty=2)
legend("topright",legend=c("Pooled","Batch 1", "Batch 2","Batch 3", "Batch 4",
  "Batch 5"), col=c("black","purple","red","blue","green","orange"),lwd=c(3,rep(1,times=5)),
  lty=c(1,1,1,2,1,2))

slopediffb1b2mod4a<- quantile(sims4a[,6] - sims4a[,7],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb1b3mod4a<- quantile(sims4a[,6] - sims4a[,8],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb1b4mod4a<- quantile(sims4a[,6] - sims4a[,9],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb1b5mod4a<- quantile(sims4a[,6] - sims4a[,10],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb2b3mod4a<- quantile(sims4a[,7] - sims4a[,8],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb2b4mod4a<- quantile(sims4a[,7] - sims4a[,9],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb2b5mod4a<- quantile(sims4a[,7] - sims4a[,10],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb3b4mod4a<- quantile(sims4a[,8] - sims4a[,9],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb3b5mod4a<- quantile(sims4a[,8] - sims4a[,10],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb4b5mod4a<- quantile(sims4a[,9] - sims4a[,10],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
## also calculate maximum posterior difference, perhaps we know more information about 
## absolute difference among batches
maxdiffmod4a<-vector()
for(i in 1:nrow(sims4a)){
  maxdiffmod4a[i]<-max(sims4a[i,6:10]) - min(sims4a[i,6:10])
}

slopediffsmod4a<-rbind(slopediffb1b2mod4a,slopediffb1b3mod4a,slopediffb1b4mod4a,slopediffb1b5mod4a,
  slopediffb2b3mod4a,slopediffb2b4mod4a,slopediffb2b5mod4a,slopediffb3b4mod4a,slopediffb3b5mod4a,
  slopediffb4b5mod4a,quantile(maxdiffmod4a,probs=c(0.025,0.975)))#,.25/2,1-0.25/2)))

rownames(slopediffsmod4a)<-c("B1-B2","B1-B3","B1-B4","B1-B5","B2-B3","B2-B4","B2-B5",
  "B3-B4","B3-B5","B4-B5","Maximum Diff")

colnames(slopediffsmod4a)<-c("2.5 %","97.5 %")

xtable(slopediffsmod4a,digits=3)
  

# b1stabmod1a<-(90-sims1a[,1])/sims1a[,6] + 8
# b2stabmod1a<-(90-sims1a[,2])/sims1a[,7] + 8
# b3stabmod1a<-(90-sims1a[,3])/sims1a[,8] + 8
# b4stabmod1a<-(90-sims1a[,4])/sims1a[,9] + 8
# b5stabmod1a<-(90-sims1a[,5])/sims1a[,10] + 8
# meanstabmod1a<-(90-sims1a[,11])/sims1a[,12] + 8
# 
# b1stabmod1b<-(90-sims1b[,1])/sims1b[,6] + 8
# b2stabmod1b<-(90-sims1b[,2])/sims1b[,7] + 8
# b3stabmod1b<-(90-sims1b[,3])/sims1b[,8] + 8
# b4stabmod1b<-(90-sims1b[,4])/sims1b[,9] + 8
# b5stabmod1b<-(90-sims1b[,5])/sims1b[,10] + 8
# meanstabmod1b<-(90-sims1b[,11])/sims1b[,12] + 8

b1stabmod4a<-(90-sims4a[,1])/sims4a[,6] + 8
b2stabmod4a<-(90-sims4a[,2])/sims4a[,7] + 8
b3stabmod4a<-(90-sims4a[,3])/sims4a[,8] + 8
b4stabmod4a<-(90-sims4a[,4])/sims4a[,9] + 8
b5stabmod4a<-(90-sims4a[,5])/sims4a[,10] + 8
meanstabmod4a<-(90-sims4a[,11])/sims4a[,12] + 8

# b1stabmod4b<-(90-sims4b[,1])/sims4b[,6] + 8
# b2stabmod4b<-(90-sims4b[,2])/sims4b[,7] + 8
# b3stabmod4b<-(90-sims4b[,3])/sims4b[,8] + 8
# b4stabmod4b<-(90-sims4b[,4])/sims4b[,9] + 8
# b5stabmod4b<-(90-sims4b[,5])/sims4b[,10] + 8
# meanstabmod4b<-(90-sims4b[,11])/sims4b[,12] + 8

chow<-c(27.5,33.5,NA,51.4,28.6,NA)

# eststabmod1a<-cbind(quantile(b1stabmod1a,prob=c(0.05)),quantile(b2stabmod1a,prob=c(0.05)),
#   quantile(b3stabmod1a,prob=c(0.05)),quantile(b4stabmod1a,prob=c(0.05)),
#   quantile(b5stabmod1a,prob=c(0.05)),quantile(meanstabmod1a,prob=c(0.05)))
# 
# eststabmod1b<-cbind(quantile(b1stabmod1b,prob=c(0.05)),quantile(b2stabmod1b,prob=c(0.05)),
#   quantile(b3stabmod1b,prob=c(0.05)),quantile(b4stabmod1b,prob=c(0.05)),
#   quantile(b5stabmod1b,prob=c(0.05)),quantile(meanstabmod1b,prob=c(0.05)))
  
eststabmod4a<-cbind(quantile(b1stabmod4a,prob=c(0.05)),quantile(b2stabmod4a,prob=c(0.05)),
  quantile(b3stabmod4a,prob=c(0.05)),quantile(b4stabmod4a,prob=c(0.05)),
  quantile(b5stabmod4a,prob=c(0.05)),quantile(meanstabmod4a,prob=c(0.05)))

# eststabmod4b<-cbind(quantile(b1stabmod4b,prob=c(0.05)),quantile(b2stabmod4b,prob=c(0.05)),
#   quantile(b3stabmod4b,prob=c(0.05)),quantile(b4stabmod4b,prob=c(0.05)),
#   quantile(b5stabmod4b,prob=c(0.05)),quantile(meanstabmod4b,prob=c(0.05)))

# eststab<-round(rbind(chow,eststabmod1a,eststabmod1b,eststabmod4a,eststabmod4b),1)
# rownames(eststab)<-c("LR","Bayes Model 1a","Bayes Model 1b","Bayes Model 4a","Bayes Model 4b")
# colnames(eststab)<-c("Batch 1","Batch 2","Batch 3","Batch 4","Batch 5","Mean")

eststab<-round(rbind(chow,eststabmod4a),1)
rownames(eststab)<-c("LR","Bayes Model 2")
colnames(eststab)<-c("Batch 1","Batch 2","Batch 3","Batch 4","Batch 5","Mean")

xtable(eststab)