setwd("/home/ed/Dropbox/MATH352/HW10")
setwd("C:\\Users\\Owner\\Dropbox\\MATH352")


### get data and CI for r
data = read.table("law_data.txt",header=F)
n = nrow(data)
Lbar = mean(data[,1])
Gbar = mean(data[,2])
r = sum((data[,1]-Lbar)*(data[,2]-Gbar))/
	sqrt( sum( (data[,1]-Lbar)^2)*sum((data[,2]-Gbar)^2))
Y = 0.5*log((1+r)/(1-r))
bound = function(Y,z,n){
	num = exp( 2* ( Y-z*sqrt(1/(n-3)))) - 1
	den = exp( 2* ( Y-z*sqrt(1/(n-3)))) + 1
	lim = num/den
	return(lim)
}
bound(Y,qnorm(0.95),n)
[1] 0.6702239
bound(Y,qnorm(0.05),n)
[1] 0.8278518

### bootstrap for percentile
B = 10000
rstar = vector("numeric",length=B)

for(i in 1:B){
	indices = sample(1:n,replace=T)
	data.star = data[indices,]
	Lbar.star = mean(data.star[,1])
	Gbar.star = mean(data.star[,2])
	rstar[i] = sum((data.star[,1]-Lbar.star)*(data.star[,2]-Gbar.star))/
		sqrt( sum( (data.star[,1]-Lbar.star)^2)*sum((data.star[,2]-Gbar.star)^2))
}

quantile(rstar,prob=c(0.05,0.95))

### BCa method

rjack = vector("numeric",length=n)

for(i in 1:n){
	data.jack = data[-i,]
	Lbar.jack = mean(data.jack[,1])
	Gbar.jack = mean(data.jack[,2])
	rjack[i] = sum((data.jack[,1]-Lbar.jack)*(data.jack[,2]-Gbar.jack))/
		sqrt( sum( (data.jack[,1]-Lbar.jack)^2)*sum((data.jack[,2]-Gbar.jack)^2))
}
	
rjack.mean = mean(rjack)

a = sum((rjack.mean - rjack)^3) / (6*(sum((rjack.mean - rjack)^2))^(3/2))

P = sum(rstar < r) / B

znaught = qnorm(P)
zalpha = qnorm(0.05)

alpha1 = pnorm(znaught + (znaught+zalpha)/(1-a*(znaught+zalpha)))
alpha2 = pnorm(znaught + (znaught-zalpha)/(1-a*(znaught-zalpha)))

quantile(rstar,prob=c(alpha1,alpha2))
