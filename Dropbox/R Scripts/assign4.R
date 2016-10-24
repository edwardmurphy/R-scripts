setwd("/home/ed/Dropbox/MATH352")

k <- c(8,9,10)
alpha <- sum(dbinom(k,10,0.5))
# this is power at null a.k.a. the alpha level

# so what happens as p increases?
alpha.star <- sum(dbinom(k,10,0.75))

# if we know underlying F, we can use this to determine actual p 
# and then determine the prob(reject null) under this condition

# get probability
F0<-2*pnorm(-1,sd=1.48)
F1<-2*pnorm(-1,sd=1.75) #prob that abs(x)>1 given sig
F2<-2*pnorm(-1,sd=2.0)
F3<-2*pnorm(-1,sd=2.5)
F4<-2*pnorm(-1,sd=3.5)

# get power (same rejection region)
B0 <- sum(dbinom(k,10,F0))
B1 <- sum(dbinom(k,10,F1))
B2 <- sum(dbinom(k,10,F2))
B3 <- sum(dbinom(k,10,F3))
B4 <- sum(dbinom(k,10,F4))

# plot power
pdf("assign4_power.pdf")
plot(c(1.48,1.75,2.0,2.5,3.5),c(B0,B1,B2,B3,B4),xlab="Sigma",
	ylab="Power",type="l",main="Power Function of Sign Test for F~N(0,Sigma^2)")
dev.off()
