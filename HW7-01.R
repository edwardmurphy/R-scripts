setwd("C:\\Users\\emurphy\\Dropbox\\MATH353")

n = 100
iter = 1000

#### Cauchy
med = vector("numeric",length=iter)
delt = vector("numeric",length=iter)

for(i in 1:iter){
	smp = rcauchy(100)
	med[i] = sort(smp)[n/2]
	num = sum( (smp-med[i])/(1+(smp-med[i])^2) )
	den = sum(  ((smp-med[i])^2-1)/((1+(smp-med[i])^2)^2)  )
	delt[i] = med[i] - num/den
}

med.norm = sqrt(n)*med
delt.norm = sqrt(n)*delt

var(med.norm)
var(delt.norm)

png("cauchy.png")
par(mfrow=c(2,1))
hist(med.norm,breaks=100,main="Histogram of normalized median for Cauchy",
	xlim=c(-6,6))
hist(delt.norm,breaks=100,main="Histogram of normalized delta for Cauchy",
	xlim=c(-6,6))
dev.off()



#### Logistic

med = vector("numeric",length=iter)
delt = vector("numeric",length=iter)

for(i in 1:iter){
	smp = rlogis(100)
	med[i] = sort(smp)[n/2]
	num = sum( (exp(smp-med[i])-1) / (exp(smp-med[i])+1))  
	den = sum( -2*exp(smp-med[i])/(exp(smp-med[i])+1)^2 )
	delt[i] = med[i] - num/den
}

med.norm = sqrt(n)*med
delt.norm = sqrt(n)*delt

var(med.norm)
var(delt.norm)

png("logistic.png")
par(mfrow=c(2,1))
hist(med.norm,breaks=100,main="Histogram of normalized median for Logistic",
	xlim=c(-6,6))
hist(delt.norm,breaks=100,main="Histogram of normalized delta for Logistic",
	xlim=c(-6,6))
dev.off()

