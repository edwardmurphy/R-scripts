#### reject rate is NOT ~binom
#### counts are ~ binom

#### appears to be a time effect? rates getting higher over time

#### need to consider prediction or tolerance intervals here
#### NOT confidence intervals

#### bootstrap the slope parameters for the size effect
#### test significance that way

#### bootstrap the CI for reject rates
#### use upper CI to sample and get pred/tol interval

############################################################################
### Initialization
############################################################################
library(tolerance)

setwd("/home/ed/Dropbox/ConsultBusiness/BiotestRejectRate")
setwd("C:\\users\\emurphy\\Dropbox\\ConsultBusiness\\BiotestRejectRate")

#data <- read.csv("RR_09JUL2013.csv",header=T,sep=",")
#data <- read.table("Biotest_rejects_09JUL13.txt",header=T)
data <- read.table("Biotest_rejects_17JUL13.txt",header=T)
as.Date(data$Date)
data <- data[order(data$Date),]

#subset data at shift
data.pre = data[1:12,]
data.post = data[13:nrow(data),]

#remove highest 3 observations
data.minout = data[fmrr<=0.04,]

attach(data)

############################################################################
### Plots
############################################################################

# plot by time (sorted)
plot(fmrr,type="l",main="Plot of Foreign Matter and Overall Reject Rates By Time",
	ylim=c(0,0.15),xlab="Time Index",ylab="Reject Rate")
points(fmrr,pch=20)
abline(h=0.05)
points(orr,type="l",lty=2,col="red")
points(orr,pch=18,col="red")
abline(h=0.10,lty=2,col="red")
#abline(v=12.5)
legend("topleft",c("Foreign Matter, Observed","Current Foreign Matter Limit","Overall,Observed","Current Overall Limit"),
	lty=c(1,1,2,2),pch=c(20,26,18,26),col=c("black","black","red","red"),cex=0.7)

#change after obs 12? 
data[12:13,]
# 5 month gap in time

# t-test
t.test(data.pre$fmrr,data.post$fmrr)
t.test(data.pre$orr,data.post$orr)

#Kruskal-Wallis test
kruskal.test(list(data.pre$fmrr,data.post$fmrr))
kruskal.test(list(data.pre$orr,data.post$orr))

# size effect?
data.n.sort <- data[order(N),]
plot(data.n.sort$N,data.n.sort$fmrr,type="l",ylim=c(0,0.12),main="Plot of Foreign Matter and Overall Reject Rates By Lot Size")
points(N,fmrr,pch=20)
points(data.n.sort$N,data.n.sort$orr,type="l",lty=2,col="red")
points(N,orr,pch=18,col="red")
legend("topleft",c("Foreign Matter","Overall"),lty=c(1,2),pch=c(20,18),col=c("black","red"))

fit.fm = lm(fmrr~N) #just above 0.05
fit.o = lm(orr~N) #yes but small R2

# is size increasing over time (so time shift/size effect confounded)?
plot(N)

# first lot is potential outlier, remove and do fit
fit2.fm = lm(data[-1,]$fmrr~data[-1,]$N)
fit2.o = lm(data[-1,]$orr~data[-1,]$N)

############################################################################
### statistics 
############################################################################

# foreign matter reject rate stats (all)
fmrr.mean = mean(fmrr)
fmrr.n = mean(N) #since different n for each lot, use mean
fmrr.sd.binom = sqrt(fmrr.mean*(1-fmrr.mean)*1/fmrr.n) #binomial sd
fmrr.sd.norm = sd(fmrr) #std error

# foreign matter reject rate stats (pre, n=12)
fmrr.mean.pre = mean(data.pre$fmrr)
fmrr.n.pre = mean(data.pre$N) #since different n for each lot, use mean
fmrr.sd.binom.pre = sqrt(fmrr.mean.pre*(1-fmrr.mean.pre)*1/fmrr.n.pre) #binomial sd
fmrr.sd.norm.pre = sd(data.pre$fmrr) #std error

# foreign matter reject rate stats (post, n=14)
fmrr.mean.post = mean(data.post$fmrr)
fmrr.n.post = mean(data.post$N) #since different n for each lot, use mean
fmrr.sd.binom.post = sqrt(fmrr.mean.post*(1-fmrr.mean.post)*1/fmrr.n.post) #binomial sd
fmrr.sd.norm.post = sd(data.post$fmrr) #std error

# overall reject rate stats (all)
orr.mean = mean(orr)
orr.n = mean(N) #since different n for each lot, use mean
orr.sd.binom = sqrt(orr.mean*(1-orr.mean)*1/orr.n) #binomial sd
orr.sd.norm = sd(orr) #std error

# overall reject rate stats (pre, n=12)
orr.mean.pre = mean(data.pre$orr)
orr.n.pre = mean(data.pre$N) #since different n for each lot, use mean
orr.sd.binom.pre = sqrt(orr.mean.pre*(1-orr.mean.pre)*1/orr.n.pre) #binomial sd
orr.sd.norm.pre = sd(data.pre$orr) #std error

# overall reject rate stats (post, n=14)
orr.mean.post = mean(data.post$orr)
orr.n.post = mean(data.post$N) #since different n for each lot, use mean
orr.sd.binom.post = sqrt(orr.mean.post*(1-orr.mean.post)*1/orr.n.post) #binomial sd
orr.sd.norm.post = sd(data.post$orr) #std error

c(fmrr.mean,fmrr.mean.pre,fmrr.mean.post)
c(fmrr.sd.binom,fmrr.sd.binom.pre,fmrr.sd.binom.post)
c(fmrr.sd.norm,fmrr.sd.norm.pre,fmrr.sd.norm.post)

c(orr.mean,orr.mean.pre,orr.mean.post)
c(orr.sd.binom,orr.sd.binom.pre,orr.sd.binom.post)
c(orr.sd.norm,orr.sd.norm.pre,orr.sd.norm.post)
 
# binomial sd is much smaller than std error
# typical? simulate a dataset with fixed p
smp = rbinom(26,round(fmrr.n,0),fmrr.mean)
phat = smp/fmrr.n
sqrt(mean(phat)*(1-mean(phat))*1/fmrr.n) 	# binomial sd
sd(phat)						# std error
# not typical, should be similiar

# add in error around p by sampling from actual data
smp = NULL
for(i in 1:26){
	smp[[i]] = rbinom(1,round(fmrr.n,0),sample(fmrr,1))
}
phat = smp/fmrr.n
sqrt(mean(phat)*(1-mean(phat))*1/fmrr.n)
sd(phat)
# there it is

############################################################################
### tolerance intervals using factors from Hahn and Meeker
############################################################################

# upper one-sided tol interval factors from Hahn and Meeker
# at 95% confidence
# proportions = 90%, 95%, 99%

fac = c(2.12,2.275,2.606)

# foreign matter
fmrr.utl.binom = matrix(fmrr.mean + fac*fmrr.sd.binom,ncol=3)
colnames(fmrr.utl.binom) = c("90%","95%","99%")
fmrr.utl.binom

fmrr.utl.norm.man = matrix(fmrr.mean + fac*fmrr.sd.norm,ncol=3) 
colnames(fmrr.utl.norm.man) = c("90%","95%","99%")
fmrr.utl.norm.man

# overall
orr.utl.binom = matrix(orr.mean + fac*orr.sd.binom,ncol=3)
colnames(orr.utl.binom) = c("90%","95%","99%")
orr.utl.binom

orr.utl.norm.man = matrix(orr.mean + fac*orr.sd.norm,ncol=3) 
colnames(orr.utl.norm.man) = c("90%","95%","99%")
orr.utl.norm.man

############################################################################
### tolerance intervals, tolerance package, all data
############################################################################

# non-parametric
fmrr.utl.np = nptol.int(fmrr,P=c(0.95,0.99),method="HM")
orr.utl.np = nptol.int(orr,P=c(0.95,0.99),method="HM")

# normal
fmrr.utl.norm = normtol.int(fmrr,P=c(0.90,0.95,0.99),method="HE") #all methods same
orr.utl.norm = normtol.int(orr,P=c(0.90,0.95,0.99),method="HE") #all methods same

fmrr.utl.norm.vec = round(t(fmrr.utl.norm[,5]),4)
orr.utl.norm.vec = round(t(orr.utl.norm[,5]),4)

utl.norm = rbind(fmrr.utl.norm.vec,orr.utl.norm.vec)
rownames(utl.norm) = c("Foreign Matter Reject Rate","Overall Reject Rate")
colnames(utl.norm) = c("90%","95%","99%")
write.csv(utl.norm,file="upper_tol_table.csv")

# plot
plottol(fmrr.utl.np,fmrr)
plottol(orr.utl.np,orr)

# normal, using data.minout
fmrr.utl.minout = normtol.int(data.minout$fmrr,P=c(0.90,0.95,0.99),method="HE") #all methods same
orr.utl.minout = normtol.int(data.minout$orr,P=c(0.90,0.95,0.99),method="HE") #all methods same

fmrr.utl.minout.vec = round(t(fmrr.utl.minout[,5]),4)
orr.utl.minout.vec = round(t(orr.utl.minout[,5]),4)

utl.minout = rbind(fmrr.utl.minout.vec,orr.utl.minout.vec)
rownames(utl.minout) = c("Foreign Matter Reject Rate","Overall Reject Rate")
colnames(utl.minout) = c("90%","95%","99%")
write.csv(utl.minout,file="upper_tol_table_remove_top_3.csv")


############################################################################
### tolerance intervals, tolerance package, subset data
############################################################################

# non-parametric
fmrr.utl.np.pre = nptol.int(data.pre$fmrr,P=c(0.90,0.95,0.99),method="HM")
orr.utl.np.pre = nptol.int(data.pre$orr,P=c(0.90,0.95,0.99),method="HM")

fmrr.utl.np.post = nptol.int(data.post$fmrr,P=c(0.90,0.95,0.99),method="HM")
orr.utl.np.post = nptol.int(data.post$orr,P=c(0.90,0.95,0.99),method="HM")

# normal
fmrr.utl.norm.pre = normtol.int(data.pre$fmrr,P=c(0.90,0.95,0.99),method="HE")
orr.utl.norm.pre = normtol.int(data.pre$orr,P=c(0.90,0.95,0.99),method="HE")

fmrr.utl.norm.post = normtol.int(data.post$fmrr,P=c(0.90,0.95,0.99),method="HE")
orr.utl.norm.post = normtol.int(data.post$orr,P=c(0.90,0.95,0.99),method="HE")

############################################################################
### SIMULATION
############################################################################

# so, if process average using only last 14 observations is the actual process mean,
# what would be the expected tolerance limit?
# assume 40% process variability

ctr = mean(data.post$fmrr)
utl.sim = as.vector(NULL)
utl = as.vector(NULL)
se = NULL
n = mean(data.post$N)
size = length(data.post$fmrr)
for(i in 1:10000){
	#simulate a dataset of number of rejected vials for 14 lots
	simdata = rbinom(14,round(n,0),runif(1,ctr-0.40*ctr, ctr+0.40*ctr))
	#get phat for each lot
	fmrr.sim = simdata/n
	#get overall phat
	p = mean(fmrr.sim)
	#get upper tolerance limit given data using tolerance package
	utl.sim = rbind(utl.sim,as.vector(normtol.int(fmrr.sim,P=c(0.90,0.95,0.99),method="HE")[,5]))
	#get observed std. error 
	se[[i]] = sqrt(p*(1-p)*1/n)
	#simulate normally distributed r.v.'s with mean = overall phat, sd = observed std. erro
	rr.sim = rnorm(10000,p,se[[i]])
	#get actual 95% and 99% quantiles
	utl = rbind(utl,as.vector(quantile(rr.sim,prob=c(0.90,0.95,0.99))))
}

apply(utl,2,mean)
apply(utl.sim,2,mean)
# 4% for fmrr

ctr = mean(data.post$orr)
utl.sim = as.vector(NULL)
utl = as.vector(NULL)
se = NULL
for(i in 1:10000){
	#simulate a dataset of no. of rejected vials for 14 lots
	simdata = rbinom(14,round(n,0),runif(1,ctr-0.40*ctr, ctr+0.40*ctr))
	#get phat for each lot
	orr.sim = simdata/n
	#get overall phat
	p = mean(orr.sim)
	#get upper tolerance limit given data using tolerance package
	utl.sim = rbind(utl.sim,as.vector(normtol.int(orr.sim,P=c(0.90,0.95,0.99),method="HE")[,5]))
	#get observed std. error 
	se[[i]] = sqrt(p*(1-p)*1/n)
	#simulate normally distributed r.v.'s with mean = overall phat, sd = observed std. erro
	rr.sim = rnorm(10000,p,se[[i]])
	#get actual 95% and 99% quantiles
	utl = rbind(utl,as.vector(quantile(rr.sim,prob=c(0.90,0.95,0.99))))
}

apply(utl,2,mean)
apply(utl.sim,2,mean)
# 8% for orr

############################################################################
### BOOTSTRAP
############################################################################

utl.sim = as.vector(NULL)
utl = as.vector(NULL)
se = NULL
n = mean(data.post$N)
size = length(data.post$fmrr)
for(i in 1:10000){
	#bootstrap
	simdata = sample(data.post$fmrr,size,replace=T)
	#get overall phat
	p = mean(simdata)
	#get upper tolerance limit given data using tolerance package
	utl.sim = rbind(utl.sim,as.vector(normtol.int(simdata,P=c(0.90,0.95,0.99),method="HE")[,5]))
	#get observed std. error 
	se[[i]] = sqrt(p*(1-p)*1/n)
	#simulate normally distributed r.v.'s with mean = overall phat, sd = observed std. erro
	rr.sim = rnorm(10000,p,se[[i]])
	#get actual 95% and 99% quantiles
	utl = rbind(utl,as.vector(quantile(rr.sim,prob=c(0.90,0.95,0.99))))
}

apply(utl,2,mean)
apply(utl.sim,2,mean)


utl.sim = as.vector(NULL)
utl = as.vector(NULL)
se = NULL
n = mean(data.post$N)
size = length(data.post$orr)
for(i in 1:10000){
	#bootstrap
	simdata = sample(data.post$orr,size,replace=T)
	#get overall phat
	p = mean(simdata)
	#get upper tolerance limit given data using tolerance package
	utl.sim = rbind(utl.sim,as.vector(normtol.int(simdata,P=c(0.90,0.95,0.99),method="HE")[,5]))
	#get observed std. error 
	se[[i]] = sqrt(p*(1-p)*1/n)
	#simulate normally distributed r.v.'s with mean = overall phat, sd = observed std. erro
	rr.sim = rnorm(10000,p,se[[i]])
	#get actual 95% and 99% quantiles
	utl = rbind(utl,as.vector(quantile(rr.sim,prob=c(0.90,0.95,0.99))))
}

apply(utl,2,mean)
apply(utl.sim,2,mean)

############################################################################
### MORE BOOTSTRAP
############################################################################

utl.sim = as.vector(NULL)
utl = as.vector(NULL)
utl.up = as.vector(NULL)
se = NULL
n = mean(data.post$N)
size = length(data.post$fmrr)
for(i in 1:10000){
	#bootstrap
	simdata = sample(data.post$fmrr,size,replace=T)
	#get overall phat
	p = mean(simdata)
	#get upper tolerance limit given data using tolerance package
	utl.sim = rbind(utl.sim,as.vector(normtol.int(simdata,P=c(0.90,0.95,0.99),method="HE")[,5]))
	#get observed std. error 
	se[[i]] = sqrt(p*(1-p)*1/n)
	#simulate normally distributed r.v.'s with mean = overall phat, sd = observed std. erro
	rr.sim = rnorm(10000,p,se[[i]])
	#get actual 95% and 99% quantiles
	utl = rbind(utl,as.vector(quantile(rr.sim,prob=c(0.90,0.95,0.99))))
	#get upper 95% CI limit for p
	p.up = p + qt(0.95,size-1)*se[[i]]
	#simulate normally distributed r.v.'s with mean = upper 95 phat, sd = observed std. erro
	rr.sim.up = rnorm(10000,p.up,se[[i]])
	#get actual 95% and 99% quantiles
	utl.up = rbind(utl.up,as.vector(quantile(rr.sim.up,prob=c(0.90,0.95,0.99))))
}

apply(utl,2,mean)
apply(utl.up,2,mean)
apply(utl.sim,2,mean)

############################################################################
### non-parametric
############################################################################

# quantiles
quants = quantile(fmrr,prob=c(0.05,0.95))
quants

#bootstrap
quant = NULL
avg = NULL
iter = 10000

for(i in 1:iter){
	rate = sample(fmrr,fmrr.n,replace=T)
	quant = rbind(quant,quantile(rate,prob=c(0.95,0.99)))
	avg[i] = mean(rate)
}

lims = apply(quant,2,mean)
bootmean = mean(avg)

p95failure = sum(fmrr > lims[[1]])/fmrr.n
p99failure = sum(fmrr > lims[[2]])/fmrr.n

p95failure
p99failure

rr = NULL
for(i in 1:iter){
	lotsize = sample(N,1)
	def = rbinom(1,lotsize,bootmean)
	rr[i] = def/lotsize
} 


##############################################################
########## CI code
##############################################################

### binomial CIs for proportion

# normal approx using indiv. lots
z = qnorm(0.95)
fm.phat = fmrr
fm.qhat = 1-fmrr

fm.ul = fmrr + z*sqrt(fm.phat*fm.qhat*1/N)

# overall
fm.phat.over = sum(fmd)/sum(N)

fm.ul.over = fm.phat.over + z*sqrt(fm.phat.over*(1-fm.phat.over)*1/sum(N))

# binomial limits for Foreign Matter
fmrr.n = mean(N)
fmrr.mean = mean(fmrr)
fmrr.sd = sqrt(fmrr.mean*(1-fmrr.mean)/fmrr.n)
fmrr.UL = fmrr.mean + 3*fmrr.sd
fmrr.LL = fmrr.mean - 3*fmrr.sd
if (fmrr.LL < 0) fmrr.LL = 0

# binomial limits for Overall
orr.n = mean(N)
orr.mean = mean(orr)
orr.sd = sqrt(orr.mean*(1-orr.mean)/orr.n)
orr.UL = orr.mean + 3*orr.sd
orr.LL = orr.mean - 3*orr.sd
if (orr.LL < 0) orr.LL = 0

# table
binlim = rbind(cbind(fmrr.LL,fmrr.UL),cbind(orr.LL,orr.UL))
colnames(binlim) = c("Lower Limit","Upper Limit")
rownames(binlim) = c("Foreign Matter","Overall")
binlim

fm.rand = rbinom(10000,round(fmrr.n,0),fmrr.UL)
p.rand = fm.rand / fmrr.n
q = quantile(p.rand,c(0.95,0.99))






