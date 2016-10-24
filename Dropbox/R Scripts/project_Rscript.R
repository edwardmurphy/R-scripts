setwd("C://Users//Owner//Documents//STAT 676")

library(xtable)
library(arm)
library(BRugs)
library(boa)

Sweave("presentation_Sweave.nw")

### Data

## this is the transpose of the data in WinBUGS
potency = structure(
.Data =   c(104.8, 102.5, 101.5, 102.4, 99.4, 96.5,
				 103.9, 101.9, 103.2, 99.6, 100.2, 98.8, 
				 103.5, 102.1, 101.9, 100.3, 99.2, 101.0, 
				 101.5, 100.3, 101.1, 100.6, 100.7, 98.4, 
				 106.1, 104.3, 101.5, 101.1, 99.4, 98.2),
			.Dim = c(6,5))
      
colnames(potency)<-c("Batch 1","Batch 2","Batch 3","Batch 4","Batch 5")
rownames(potency)<-c("0","3","6","9","12","18")

potencyvec<-as.vector(potency)
t<-rep(c(0,3,6,9,12,18),times=5)
#center t so beta0 and beta1 independent
tstar<-t-mean(t)
batch<-as.factor(rep(1:5,each=6))

## Create dataframe with all data
stab<-cbind(potencyvec,t,tstar,batch)
stab<-data.frame(stab)
stab$batch<-as.factor(stab$batch)

## Fit pooled for initial estimates of mu.beta0 and mu.beta1
lm.overall<-lm(potencyvec~t)

## ANCOVA fit with different slopes/intercepts among batches
lm.ANCOVA<-lm(potencyvec~t+batch+t:batch,data=stab)

lm.b1<-lm(potencyvec~t,data=subset(stab,batch==1))
lm.b2<-lm(potencyvec~t,data=subset(stab,batch==2))
lm.b3<-lm(potencyvec~t,data=subset(stab,batch==3))
lm.b4<-lm(potencyvec~t,data=subset(stab,batch==4))
lm.b5<-lm(potencyvec~t,data=subset(stab,batch==5))

#coefficients for slope/intercept from lm for separate batches
indfitcoefs<-cbind(lm.b1$coef,lm.b2$coef,lm.b3$coef,lm.b4$coef,lm.b5$coef)

######################## Bayesian Models in WinBUGS

############# simulations 

#inits when tau is stochastic
inits1<-list(tau.g = 1 , tau.beta0 = 0.1, tau.beta1 = 0.1)#, mu.beta0 = 101, mu.beta1 = -0.3)
inits2<-list(tau.g = 1 , tau.beta0 = 1, tau.beta1 = 1)#, mu.beta0 = 101, mu.beta1 = -0.3)
inits3<-list(tau.g = 1 , tau.beta0 = 10, tau.beta1 = 10)#, mu.beta0 = 101, mu.beta1 = -0.3)
inits<-list(inits1,inits2,inits3)

parameters<-c("beta0[]","beta1[]","mu.beta0","mu.beta1","tau.g","tau.beta0","tau.beta1",
"sig2.g","sig2.beta0","sig2.beta1")


##### Model 1

#non-informative priors
#post means: tau.g = 1.0, tau.beta0 = 166.7, tau.beta1 = 119.2, mu.beta0 = 101.2,
#   mu.beta1 = -0.3

#burn-in is 1000, variance (and precision) are not very stable
#with 15000 iterations/chain slopes are accurate to three decimal places, ints to 
#  2 decimal places


# sim1burn<-bugs(data="project_data.txt",inits,parameters,"project_bugsmodel.txt",
#   debug=T,
#   n.chains=3,n.iter=2000,n.burnin=1000,n.thin=1,
#   bugs.seed=floor(runif(1,1,10000)))

sim1<-bugs(data="project_data.txt",inits,parameters,"project_bugsmodel.txt",
  #debug=T,
  n.chains=3,n.iter=16000,n.burnin=1000,n.thin=1,
  bugs.seed=floor(runif(1,1,10000)))

sims<-sim1$sims.matrix

##### Model 2

###uniform priors on mu.beta0 and mu.beta1 within appropriate ranges, i.e. more
## informative priors on these parameters

#burn-in is 1000, much faster convergence
#with 15000 iterations/chain slopes are accurate to three decimal places, ints to 
#  2 decimal places

# sim2burn<-bugs(data="project_data.txt",inits,parameters,"project_bugsmodel2.txt",
#   debug=T,
#   n.chains=3,n.iter=2000,n.burnin=1000,n.thin=1,
#   bugs.seed=floor(runif(1,1,10000)))

sim2<-bugs(data="project_data.txt",inits,parameters,"project_bugsmodel2.txt",
  #debug=T,
  n.chains=3,n.iter=16000,n.burnin=1000,n.thin=1,
  bugs.seed=floor(runif(1,1,10000)))

sims2<-sim2$sims.matrix

### models 1 and 2 give nearly identical results, the real difference lies in the
### choice of prior for the variance (see model 4)

##### Model 4

###uniform priors on mu.beta0 and mu.beta1 within appropriate ranges, i.e. more
## informative priors on these parameters, uniform on all standard deviation terms and 
## tau are logical

#inits when tau is stochastic
inits1<-list(sig.g = 0.1 , sig.beta0 = 0.1, sig.beta1 = 0.01)#, mu.beta0 = 101, mu.beta1 = -0.3)
inits2<-list(sig.g = 1 , sig.beta0 = 1, sig.beta1 = 0.1)#, mu.beta0 = 101, mu.beta1 = -0.3)
inits3<-list(sig.g = 10 , sig.beta0 = 5, sig.beta1 = 1)#, mu.beta0 = 101, mu.beta1 = -0.3)
inits<-list(inits1,inits2,inits3)

# sim4burn<-bugs(data="project_data.txt",inits,parameters,"project_bugsmodel4.txt",
#   debug=T,
#   n.chains=3,n.iter=2000,n.burnin=1000,n.thin=1,
#   bugs.seed=floor(runif(1,1,10000)))

sim4<-bugs(data="project_data.txt",inits,parameters,"project_bugsmodel4.txt",
  #debug=T,
  n.chains=3,n.iter=16000,n.burnin=1000,n.thin=1,
  bugs.seed=floor(runif(1,1,10000)))

sims4<-sim4$sims.matrix

## Original analysis was comparing priors on means (model 1 v 2), but perhaps it is 
## more apropos to compare priors on variance terms (e.g. model 2 v 4). Analyses below
## were again originally written to compare means, but changes may be made to compare
## variances.  thus, may not be consistent if saved in the middle of analyses

#### overlay kernel density plots of posterior sig2.beta1 (difference in slopes)

x11()
plot(density(sims[,18]),xlim=c(0,0.5),main="Posterior Density of Batch Slope Variance",
  xlab="")
lines(density(sims2[,18]),lty=2)
lines(density(sims4[,18]),lty=3)
legend("topright",legend=c("Model 1","Model 2","Model 4"),lty=1:3)

##overlay kernel density plots of posterior overall slope with individual slopes, 
##looks bi-modal, maybe tri-modal? perhaps our batches have unknown differences 
##which are leading to 2 (or 3) different degradation paths
# no differences between models though
x11()
par(mfrow=c(3,1))
#Model 1
plot(density(sims[,12]),lwd=3,xlim=c(-1,1),
  main="Pooled and Individual Density Plots for Slope, Model 1",
  xlab="Posterior Slope",ylim=c(0,7))
lines(density(sims[,6]),col="purple")
lines(density(sims[,7]),col="red")
lines(density(sims[,8]),col="blue",lty=2)
lines(density(sims[,9]),col="green")
lines(density(sims[,10]),col="orange",lty=2)
legend("topright",legend=c("Pooled","Batch 1", "Batch 2","Batch 3", "Batch 4",
  "Batch 5"), col=c("black","purple","red","blue","green","orange"),lwd=c(3,rep(1,times=5)),
  lty=c(1,1,1,2,1,2))
#Model 2
plot(density(sims2[,12]),lwd=3,xlim=c(-1,1),
  main="Pooled and Individual Density Plots for Slope, Model 2",
  xlab="Posterior Slope",ylim=c(0,7))
lines(density(sims2[,6]),col="purple")
lines(density(sims2[,7]),col="red")
lines(density(sims2[,8]),col="blue",lty=2)
lines(density(sims2[,9]),col="green")
lines(density(sims2[,10]),col="orange",lty=2)
legend("topright",legend=c("Pooled","Batch 1", "Batch 2","Batch 3", "Batch 4",
  "Batch 5"), col=c("black","purple","red","blue","green","orange"),lwd=c(3,rep(1,times=5)),
  lty=c(1,1,1,2,1,2))
#Model 4
plot(density(sims4[,12]),lwd=3,xlim=c(-1,1),
  main="Pooled and Individual Density Plots for Slope, Model 4",
  xlab="Posterior Slope",ylim=c(0,7))
lines(density(sims4[,6]),col="purple")
lines(density(sims4[,7]),col="red")
lines(density(sims4[,8]),col="blue",lty=2)
lines(density(sims4[,9]),col="green")
lines(density(sims4[,10]),col="orange",lty=2)
legend("topright",legend=c("Pooled","Batch 1", "Batch 2","Batch 3", "Batch 4",
  "Batch 5"), col=c("black","purple","red","blue","green","orange"),lwd=c(3,rep(1,times=5)),
  lty=c(1,1,1,2,1,2))

##boxplots??

#calculate quantiles of pairwise differences in posterior batch slopes
# if does not contain zero, then a significant difference...depends on what 
# sort of probability you want. unknown how related to p = 0.25 criterion
# differences noted
# model 1
slopediffb1b2mod1<- quantile(sims[,6] - sims[,7],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb1b3mod1<- quantile(sims[,6] - sims[,8],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb1b4mod1<- quantile(sims[,6] - sims[,9],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb1b5mod1<- quantile(sims[,6] - sims[,10],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb2b3mod1<- quantile(sims[,7] - sims[,8],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb2b4mod1<- quantile(sims[,7] - sims[,9],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb2b5mod1<- quantile(sims[,7] - sims[,10],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb3b4mod1<- quantile(sims[,8] - sims[,9],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb3b5mod1<- quantile(sims[,8] - sims[,10],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb4b5mod1<- quantile(sims[,9] - sims[,10],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
## also calculate maximum posterior difference, perhaps we know more information about 
## absolute difference among batches
maxdiffmod1<-vector()
for(i in 1:nrow(sims)){
  maxdiffmod1[i]<-max(sims[i,6:10]) - min(sims[i,6:10])
}

slopediffb1b2mod2<- quantile(sims2[,6] - sims2[,7],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb1b3mod2<- quantile(sims2[,6] - sims2[,8],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb1b4mod2<- quantile(sims2[,6] - sims2[,9],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb1b5mod2<- quantile(sims2[,6] - sims2[,10],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb2b3mod2<- quantile(sims2[,7] - sims2[,8],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb2b4mod2<- quantile(sims2[,7] - sims2[,9],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb2b5mod2<- quantile(sims2[,7] - sims2[,10],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb3b4mod2<- quantile(sims2[,8] - sims2[,9],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb3b5mod2<- quantile(sims2[,8] - sims2[,10],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb4b5mod2<- quantile(sims2[,9] - sims2[,10],probs=c(0.025,0.975))#,.25/2,1-0.25/2))

maxdiffmod2<-vector()
for(i in 1:nrow(sims2)){
  maxdiffmod2[i]<-max(sims2[i,6:10]) - min(sims2[i,6:10])
}

slopediffb1b2mod4<- quantile(sims4[,6] - sims4[,7],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb1b3mod4<- quantile(sims4[,6] - sims4[,8],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb1b4mod4<- quantile(sims4[,6] - sims4[,9],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb1b5mod4<- quantile(sims4[,6] - sims4[,10],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb2b3mod4<- quantile(sims4[,7] - sims4[,8],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb2b4mod4<- quantile(sims4[,7] - sims4[,9],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb2b5mod4<- quantile(sims4[,7] - sims4[,10],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb3b4mod4<- quantile(sims4[,8] - sims4[,9],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb3b5mod4<- quantile(sims4[,8] - sims4[,10],probs=c(0.025,0.975))#,.25/2,1-0.25/2))
slopediffb4b5mod4<- quantile(sims4[,9] - sims4[,10],probs=c(0.025,0.975))#,.25/2,1-0.25/2))

maxdiffmod4<-vector()
for(i in 1:nrow(sims4)){
  maxdiffmod4[i]<-max(sims4[i,6:10]) - min(sims4[i,6:10])
}

slopediffsmod1<-rbind(slopediffb1b2mod1,slopediffb1b3mod1,slopediffb1b4mod1,slopediffb1b5mod1,
  slopediffb2b3mod1,slopediffb2b4mod1,slopediffb2b5mod1,slopediffb3b4mod1,slopediffb3b5mod1,
  slopediffb4b5mod1,quantile(maxdiffmod1,probs=c(0.025,0.975)))#,.25/2,1-0.25/2)))
  
slopediffsmod2<-rbind(slopediffb1b2mod2,slopediffb1b3mod2,slopediffb1b4mod2,slopediffb1b5mod2,
  slopediffb2b3mod2,slopediffb2b4mod2,slopediffb2b5mod2,slopediffb3b4mod2,slopediffb3b5mod2,
  slopediffb4b5mod2,quantile(maxdiffmod2,probs=c(0.025,0.975)))#,.25/2,1-0.25/2)))

slopediffsmod4<-rbind(slopediffb1b2mod4,slopediffb1b3mod4,slopediffb1b4mod4,slopediffb1b5mod4,
  slopediffb2b3mod4,slopediffb2b4mod4,slopediffb2b5mod4,slopediffb3b4mod4,slopediffb3b5mod4,
  slopediffb4b5mod4,quantile(maxdiffmod4,probs=c(0.025,0.975)))#,.25/2,1-0.25/2)))

slopediffs<-cbind(slopediffsmod1,slopediffsmod2,slopediffsmod4)

rownames(slopediffs)<-c("B1-B2","B1-B3","B1-B4","B1-B5","B2-B3","B2-B4","B2-B5",
  "B3-B4","B3-B5","B4-B5","Maximum Diff")

colnames(slopediffs)<-c("2.5%, Model 1","97.5%, Model 1", "2.5%, Model 2", "97.5%, Model 2",
  "2.5%, Model 4", "97.5%, Model 4")

##calculate posterior estimates of mean shelf-life (expectations)
## not predictive (no new draws)
b1stabmod1<-(90-sims[,1])/sims[,6] + 8
b2stabmod1<-(90-sims[,2])/sims[,7] + 8
b3stabmod1<-(90-sims[,3])/sims[,8] + 8
b4stabmod1<-(90-sims[,4])/sims[,9] + 8
b5stabmod1<-(90-sims[,5])/sims[,10] + 8
meanstabmod1<-(90-sims[,11])/sims[,12] + 8

b1stabmod2<-(90-sims2[,1])/sims2[,6] + 8
b2stabmod2<-(90-sims2[,2])/sims2[,7] + 8
b3stabmod2<-(90-sims2[,3])/sims2[,8] + 8
b4stabmod2<-(90-sims2[,4])/sims2[,9] + 8
b5stabmod2<-(90-sims2[,5])/sims2[,10] + 8
meanstabmod2<-(90-sims2[,11])/sims2[,12] + 8

b1stabmod4<-(90-sims4[,1])/sims4[,6] + 8
b2stabmod4<-(90-sims4[,2])/sims4[,7] + 8
b3stabmod4<-(90-sims4[,3])/sims4[,8] + 8
b4stabmod4<-(90-sims4[,4])/sims4[,9] + 8
b5stabmod4<-(90-sims4[,5])/sims4[,10] + 8
meanstabmod4<-(90-sims4[,11])/sims4[,12] + 8

##5% lower credible limit shows shrinkage...some higher, some lower. cool!!
chow<-c(27.5,33.5,NA,51.4,28.6)

eststabmod1<-cbind(quantile(b1stabmod1,prob=c(0.05)),quantile(b2stabmod1,prob=c(0.05)),
  quantile(b3stabmod1,prob=c(0.05)),quantile(b4stabmod1,prob=c(0.05)),
  quantile(b5stabmod1,prob=c(0.05)))

eststabmod2<-cbind(quantile(b1stabmod2,prob=c(0.05)),quantile(b2stabmod2,prob=c(0.05)),
  quantile(b3stabmod2,prob=c(0.05)),quantile(b4stabmod2,prob=c(0.05)),
  quantile(b5stabmod2,prob=c(0.05)))

eststabmod4<-cbind(quantile(b1stabmod4,prob=c(0.05)),quantile(b2stabmod4,prob=c(0.05)),
  quantile(b3stabmod4,prob=c(0.05)),quantile(b4stabmod4,prob=c(0.05)),
  quantile(b5stabmod4,prob=c(0.05)))

eststab<-round(rbind(chow,eststabmod1,eststabmod2,eststabmod4),1)
rownames(eststab)<-c("LR","Bayes Model 1","Bayes Model 2","Bayes Model 4")
colnames(eststab)<-c("Batch 1","Batch 2","Batch 3","Batch 4","Batch 5","Mean")

xtable(eststab)

#plot all posterior regression lines with overall mean
# note curvature and skewness, so this is why credible set is low
x11()
plot(1,2,type="n",xlim=c(-8,62),ylim=c(85,110),xlab="Time (months)",ylab="Potency (%)",xaxs="i",xaxt="n")
for (i in 1:nrow(sims)){
  abline(a = sims[i,1], b = sims[i,6], col="grey")
}
abline(h = 90, lty=2)
abline(a = mean(sims[,1]), b = mean(sims[,6]), lwd=3)

axis(1,at=c(-8,2,12,22,32,42,52,62),labels=c(0,10,20,30,40,50,60,70))

#because of skewness, highest posterior density region more appropriate??
hdrs<-c(boa.hpd(b1stab,0.95)[1],boa.hpd(b2stab,0.95)[1],boa.hpd(b3stab,0.95)[1],
  boa.hpd(b4stab,0.95)[1],boa.hpd(b5stab,0.95)[1])

eststab<-rbind(hdrs,eststab)



