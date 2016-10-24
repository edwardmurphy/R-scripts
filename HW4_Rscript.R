setwd("C://Users//Owner//Documents//STAT 676")
Sweave("HW4_Sweave.nw")

library(xtable)
library(arm)
library(BRugs)

runningmeanplots<-function(data,zoom){
  par(mfrow=c(length(data)-1,1))
  if (zoom)
  for (i in 1:(length(data)-1)){
    plot(cumsum(data[[i]])/seq(along=data[[i]]),type="l",ylab="Running Mean",main=names(data)[[i]],ylim=c(0.99*(mean(data[[i]])),1.01*mean(data[[i]])))
    abline(h=mean(data[[i]]))
  }
  else
  for (i in 1:(length(data)-1)){
    plot(cumsum(data[[i]])/seq(along=data[[i]]),type="l",ylab="Running Mean",main=names(data)[[i]])
    abline(h=mean(data[[i]]))
  }
}

halfwidths<-function(data){
  hw<-rep(0,times=length(data))
  for(i in 1:length(data)){
    postdata<-data[[i]]
    m=length(postdata)
    n=floor(m^0.5)
    k=floor(m/n)
    gmean=mean(postdata)
    bmean=rep(0,times=k)
    for (j in 1:k){
       bmean[j]=(1/n)*sum(postdata[((j-1)*n+1):(j*n)])
    }
    varhat=(1/(k*(k-1)))*sum((bmean-gmean)^2)
    se=sqrt(varhat)
    hw[i]<-qt(0.975,k-1)*se
  }
  return(hw)
}

###Problem 1a
#select n.iter the same for all three subproblems
n.iter=36000
# selected so that halfwidth of mu, theta <= 0.004, so estimate is accurate to 2 decimal places

y<-c(1.2,2.4,1.3,1.3,0.0,1.0,1.8,0.8,4.6,1.4)
n<-length(y)

data<-list("n","y")
inits<-function(){
  list(mu=0,tau=1)
}
parameters<-c("mu","tau","x")
P1a.sim<-bugs(data,inits,parameters,"HW4_1a_model.txt",
  #debug=T,
  n.chains=1,n.iter=n.iter,n.burnin=1000,n.thin=1,bugs.seed=floor(runif(1,1,10000)))

P1a.post<-P1a.sim$sims.list
P1a.sum<-cbind(P1a.sim$summary,halfwidths(P1a.post))
colnames(P1a.sum)<-c(colnames(P1a.sim$summary),"MC error")

runningmeanplots(P1a.post,zoom=T)

#use x for post prob that diff is postive, or alternatively
#     sum(P1a.post$mu>0)/length(P1a.post$mu)

###Problem 1b
P1b.sim<-bugs(data,inits,parameters,"HW4_1b_model.txt",
  #debug=T,
  n.chains=1,n.iter=n.iter,n.burnin=1000,n.thin=1,bugs.seed=floor(runif(1,1,10000)))

P1b.post<-P1b.sim$sims.list
P1b.sum<-cbind(P1b.sim$summary,halfwidths(P1b.post))
colnames(P1b.sum)<-c(colnames(P1b.sim$summary),"MC error")

runningmeanplots(P1b.post,zoom=T)
## mu converges faster than normal like, but tau is much slower to converge

###Problem 1c
y<-rep.int(1,times=9) #all differences positive after removing 0.0
n<-length(y)

data<-list("n","y")
inits<-function(){
  theta <- rbeta(1,1,1)
}
parameters<-c("theta","x")
P1c.sim<-bugs(data,inits,parameters,"HW4_1c_model.txt",
  #debug=T,
  n.chains=1,n.iter=n.iter,n.burnin=1000,n.thin=1,bugs.seed=floor(runif(1,1,10000)))

P1c.post<-P1c.sim$sims.list
P1c.sum<-cbind(P1c.sim$summary,halfwidths(P1c.post))
colnames(P1c.sum)<-c(colnames(P1c.sim$summary),"MC error")

runningmeanplots(P1c.post,zoom=T)

## Note big difference in tau estimates b/w parts a & b

postprob<-c(P1a.sum[3,1],P1b.sum[3,1],P1c.sum[2,1])
paramest<-round(c(P1a.sum[1,1],P1b.sum[1,1],P1c.sum[1,1]),2)
paramcredlow<-round(c(P1a.sum[1,3],P1b.sum[1,3],P1c.sum[1,3]),2)
paramcredhigh<-round(c(P1a.sum[1,7],P1b.sum[1,7],P1c.sum[1,7]),2)

p1tab<-cbind(postprob,paramest,paramcredlow,paramcredhigh)
rownames(p1tab)<-c("Normal","Student's t","Binomial")
  
tauest<-c(P1a.sum[2,1],P1b.sum[2,1])
taucredlow<-c(P1a.sum[2,3],P1b.sum[2,3])
taucredhigh<-c(P1a.sum[2,7],P1b.sum[2,7])

tautab<-cbind(tauest,taucredlow,taucredhigh)
rownames(tautab)<-c("Normal","Student's t")




###Problem 2
dunnock<-c(22.0,23.9,20.9,23.8,25.0,24.0,21.7,23.8,22.8,23.1,23.5,23.0,23.1,23.0)
warbler<-c(23.2,22.0,22.2,21.2,21.6,21.6,21.9,22.0,22.9,22.8)

boxplot(dunnock,warbler,names=c("Dunnock","Reed Warbler"))
qqnorm(scale(dunnock),main="Normal Q-Q Plot for Dunnock")
abline(0,1)
qqnorm(scale(warbler),main="Normal Q-Q Plot for Reed Warbler")
abline(0,1)

shapiro.test(dunnock)
shapiro.test(warbler)

t.test(dunnock,warbler,var.equal=T)
t.test(dunnock,warbler)

tpooledint<-c(t.test(dunnock,warbler,var.equal=T)$conf.int[1],t.test(dunnock,warbler,var.equal=T)$conf.int[2])
tsepint<-c(t.test(dunnock,warbler)$conf.int[1],t.test(dunnock,warbler)$conf.int[2])
tints<-rbind(tpooledint,tsepint)


##Data names in WinBUGS model
xd<-dunnock
nd<-length(xd)
xr<-warbler
nr<-length(xr)

n.iter=36000

##Bayes model with pooled variance
data<-list("nd","nr","xd","xr")
inits<-function(){
  list(mud=22.2,mur=22.2,tau=1)
}
parameters<-c("mud","mur","tau","mudelta")
P2a.sim<-bugs(data,inits,parameters,"HW4_2a_model.txt",
  #debug=T,
  n.chains=1,n.iter=n.iter,n.burnin=1000,n.thin=1,bugs.seed=floor(runif(1,1,10000)))
  
  
P2a.post<-P2a.sim$sims.list
P2a.sum<-cbind(P2a.sim$summary,halfwidths(P2a.post))
colnames(P2a.sum)<-c(colnames(P2a.sim$summary),"MC error")
P2a.sum

##Bayes model with unpooled variance
inits<-function(){
list(mud=22.2,mur=22.2,taud=1,taur=1)
}
parameters<-c("mud","mur","taud","taur","mudelta")
P2b.sim<-bugs(data,inits,parameters,"HW4_2b_model.txt",
  #debug=T,
  n.chains=1,n.iter=n.iter,n.burnin=1000,n.thin=1,bugs.seed=floor(runif(1,1,10000)))
P2b.post<-P2b.sim$sims.list
P2b.sum<-cbind(P2b.sim$summary,halfwidths(P2b.post))
colnames(P2b.sum)<-c(colnames(P2b.sim$summary),"MC error")
P2b.sum

poolcred<-c(P2a.sum[4,3],P2a.sum[4,7])
sepcred<-c(P2b.sum[5,3],P2b.sum[5,7])
creds<-rbind(poolcred,sepcred)

##informative priors
n.iter=26000
## requires less iterations than informative prior


##Bayes model with pooled variance, informative normal prior for mean
inits<-function(){
  list(mud=22.2,mur=22.2,tau=1)
}
parameters<-c("mud","mur","tau","mudelta")
P2c.sim<-bugs(data,inits,parameters,"HW4_2c_model.txt",
  #debug=T,
  n.chains=1,n.iter=n.iter,n.burnin=1000,n.thin=1,bugs.seed=floor(runif(1,1,10000)))
P2c.post<-P2c.sim$sims.list
P2c.sum<-cbind(P2c.sim$summary,halfwidths(P2c.post))
colnames(P2c.sum)<-c(colnames(P2c.sim$summary),"MC error")
P2c.sum

##Bayes model with unpooled variance, informative normal prior for mean
inits<-function(){
list(mud=22.2,mur=22.2,taud=1,taur=1)
}
parameters<-c("mud","mur","taud","taur","mudelta")
P2d.sim<-bugs(data,inits,parameters,"HW4_2d_model.txt",
  #debug=T,
  n.chains=1,n.iter=n.iter,n.burnin=1000,n.thin=1,bugs.seed=floor(runif(1,1,10000)))
P2d.post<-P2d.sim$sims.list
P2d.sum<-cbind(P2d.sim$summary,halfwidths(P2d.post))
colnames(P2d.sum)<-c(colnames(P2d.sim$summary),"MC error")
P2d.sum

informpoolcred<-c(P2c.sum[4,3],P2c.sum[4,7])
informsepcred<-c(P2d.sum[5,3],P2d.sum[5,7])
informcreds<-rbind(informpoolcred,informsepcred)


