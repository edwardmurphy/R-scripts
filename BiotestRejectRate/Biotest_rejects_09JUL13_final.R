############################################################################
### Initialization
############################################################################
library(tolerance)

setwd("/home/ed/Dropbox/ConsultBusiness/BiotestRejectRate")
setwd("C:\\users\\emurphy\\Dropbox\\ConsultBusiness\\BiotestRejectRate")

data <- read.table("Biotest_rejects_09JUL13.txt",header=T)
as.Date(data$Date)
data <- data[order(data$Date),]

#subset data at shift
data.pre = data[1:12,]
data.post = data[13:nrow(data),]

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
legend("topleft",c("Foreign Matter, Observed","Current Foreign Matter Limit","Overall,Observed","Current Overall Limit"),
	lty=c(1,1,2,2),pch=c(20,26,18,26),col=c("black","black","red","red"),cex=0.7)

# plot by time (sorted), add proposed limits of 4% and 9%
plot(fmrr,type="l",main="Plot of Foreign Matter and Overall Reject Rates By Time",
	ylim=c(0,0.15),xlab="Time Index",ylab="Reject Rate")
points(fmrr,pch=20)
abline(h=0.04)
points(orr,type="l",lty=2,col="red")
points(orr,pch=18,col="red")
abline(h=0.09,lty=2,col="red")
legend("topleft",c("Foreign Matter, Observed","Proposed Foreign Matter Limit","Overall,Observed","Proposed Overall Limit"),
	lty=c(1,1,2,2),pch=c(20,26,18,26),col=c("black","black","red","red"),cex=0.7)



# t-test
t.test(data.pre$fmrr,data.post$fmrr)
t.test(data.pre$orr,data.post$orr)

# size effect?
data.n.sort <- data[order(N),]
plot(data.n.sort$N,data.n.sort$fmrr,type="l",ylim=c(0,0.12),main="Plot of Foreign Matter and Overall Reject Rates By Lot Size")
points(N,fmrr,pch=20)
points(data.n.sort$N,data.n.sort$orr,type="l",lty=2,col="red")
points(N,orr,pch=18,col="red")
legend("topleft",c("Foreign Matter","Overall"),lty=c(1,2),pch=c(20,18),col=c("black","red"))

fit.fm = lm(fmrr~N) #just above 0.05
fit.o = lm(orr~N) #yes but small R2

############################################################################
### statistics, manually calculated tolerance limits
############################################################################

# foreign matter reject rate stats
fmrr.mean = mean(fmrr)
fmrr.n = mean(N) #since different n for each lot, use mean
fmrr.sd.binom = sqrt(fmrr.mean*(1-fmrr.mean)*1/fmrr.n) #binomial sd
fmrr.sd.norm = sd(fmrr) #std error

# overall reject rate stats
orr.mean = mean(orr)
orr.n = mean(N) #since different n for each lot, use mean
orr.sd.binom = sqrt(orr.mean*(1-orr.mean)*1/orr.n) #binomial sd
orr.sd.norm = sd(orr) #std error

############################################################################
### tolerance intervals, tolerance package, all data
############################################################################

fmrr.utl.norm = normtol.int(fmrr,P=c(0.90,0.95,0.99),method="HE")
orr.utl.norm = normtol.int(orr,P=c(0.90,0.95,0.99),method="HE")

fmrr.utl.norm.vec = round(t(fmrr.utl.norm[,5]),4)
orr.utl.norm.vec = round(t(orr.utl.norm[,5]),4)

############################################################################
### tolerance intervals, tolerance package, subset data
############################################################################

fmrr.utl.norm.pre = normtol.int(data.pre$fmrr,P=c(0.90,0.95,0.99),method="HE")
orr.utl.norm.pre = normtol.int(data.pre$orr,P=c(0.90,0.95,0.99),method="HE")

fmrr.utl.norm.post = normtol.int(data.post$fmrr,P=c(0.90,0.95,0.99),method="HE")
orr.utl.norm.post = normtol.int(data.post$orr,P=c(0.90,0.95,0.99),method="HE")





