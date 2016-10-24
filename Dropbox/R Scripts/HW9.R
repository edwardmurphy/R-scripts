library(xtable)

#laptop
setwd("C:\\Users\\Owner\\Dropbox\\MATH352\\HW9")
#Ubuntu
setwd("/home/ed/Dropbox/MATH352/HW9")

data <- read.table("HW9.txt",header=F)

### model matrix
lot1 = c(rep(1,9),rep(0,18))
lot2 = c(rep(0,9),rep(1,9),rep(0,9))
lot3 = c(rep(0,18),rep(1,9))
hours = data[,2]

X = cbind(lot1,lot2,lot3,hours)
Y = data[,3]

n = length(Y)
p = ncol(X)

### contrasts

AvB = rbind(1,-1,0,0)
AvC = rbind(1,0,-1,0)
BvC = rbind(0,1,-1,0)

## matrices which will be used frequently
D = solve(t(X) %*% X) 				#(X'X)^-1
hat.mat = D%*%t(X)				#(X'X)^-1 * X'
proj = X %*% hat.mat				# X * (X'X)^-1 * X'

### linear fit and CIs

beta.ls = solve(t(X)%*%X) %*% t(X) %*% Y
sigsq.ls = 1/(n-p) * t(Y) %*% (diag(1,nrow(proj)) - proj) %*% Y
# compare to summary(lm(Y~X-1))



ci = function(contrast,betahat, lowercp, uppercp, var, tXX.inv){
	low = t(contrast)%*%betahat - lowercp*sqrt(var*t(contrast)%*%tXX.inv%*%contrast)
	up = t(contrast)%*%betahat - uppercp*sqrt(var*t(contrast)%*%tXX.inv%*%contrast)
	return(rbind(low,up))
}

avb.t = ci(AvB,beta.ls,qt(0.975,n-p),qt(0.025,n-p),sigsq.ls,D)
avc.t = ci(AvC,beta.ls,qt(0.975,n-p),qt(0.025,n-p),sigsq.ls,D)
bvc.t = ci(BvC,beta.ls,qt(0.975,n-p),qt(0.025,n-p),sigsq.ls,D)

ci.t = cbind(avb.t,avc.t,bvc.t)
colnames(ci.t) = c("A vs. B","A vs. C","B vs. C")
rownames(ci.t) = c("Lower Bound","Upper Bound")

### bootstrap residuals
r.hat = Y - X%*%beta.ls
e.hat = r.hat - mean(r.hat)

B = 10000
XB = X%*%beta.ls
avb.t.star = vector("numeric",B)
avc.t.star = vector("numeric",B)
bvc.t.star = vector("numeric",B)

for(i in 1:B){
	e.star = sample(e.hat,n,replace=T)
	y.star = XB + e.star
	beta.star = hat.mat%*%y.star
	sigsq.star = 1/(n-p)*t(y.star - X%*%beta.star)%*%(y.star - X%*%beta.star)
	
	avb.t.star[i] = t(AvB)%*%(beta.star-beta.ls)/sqrt(sigsq.star*t(AvB)%*% D %*% AvB)
	avc.t.star[i] = t(AvC)%*%(beta.star-beta.ls)/sqrt(sigsq.star*t(AvC)%*% D %*% AvC)
	bvc.t.star[i] = t(BvC)%*%(beta.star-beta.ls)/sqrt(sigsq.star*t(BvC)%*% D %*% BvC)
}

avb.res = ci(AvB,beta.ls,quantile(avb.t.star,prob=0.975),
	quantile(avb.t.star,prob=0.025),sigsq.ls,D)
avc.res = ci(AvC,beta.ls,quantile(avc.t.star,prob=0.975),
	quantile(avc.t.star,prob=0.025),sigsq.ls,D)
bvc.res = ci(BvC,beta.ls,quantile(bvc.t.star,prob=0.975),
	quantile(bvc.t.star,prob=0.025),sigsq.ls,D)

ci.res = cbind(avb.res,avc.res,bvc.res)
colnames(ci.res) = c("A vs. B","A vs. C","B vs. C")
rownames(ci.res) = c("Lower Bound","Upper Bound")

### bootstrap pairs
avb.t.star2 = vector("numeric",B)
avc.t.star2 = vector("numeric",B)
bvc.t.star2 = vector("numeric",B)

for(i in 1:B){
	indices = sample(1:n,replace=T)
	data.star = data[indices,]
	X.star = matrix(0,nrow=n,ncol=p)
	X.star[,p] = data.star[,2]
	for(j in 1:n){
		if(data.star[j,1]==1) 
			X.star[j,1] = 1
		else if(data.star[j,1]==2) 
			X.star[j,2] = 1
		else if(data.star[j,1]==3)
			X.star[j,3] = 1
	}
	D.star = solve(t(X.star) %*% X.star)
	hat.mat.star = D.star%*%t(X.star)				
	proj.star = X.star %*% hat.mat.star

	y.star = data.star[,3]
	
	beta.star = hat.mat.star %*% y.star
	sigsq.star = 1/(n-p)*t(y.star)%*%(diag(1,n) - proj.star)%*%y.star

	avb.t.star2[i] = t(AvB)%*%(beta.star-beta.ls)/sqrt(sigsq.star*t(AvB)%*% D.star %*% AvB)
	avc.t.star2[i] = t(AvC)%*%(beta.star-beta.ls)/sqrt(sigsq.star*t(AvC)%*% D.star %*% AvC)
	bvc.t.star2[i] = t(BvC)%*%(beta.star-beta.ls)/sqrt(sigsq.star*t(BvC)%*% D.star %*% BvC)
}

avb.pair = ci(AvB,beta.ls,quantile(avb.t.star2,prob=0.975),
	quantile(avb.t.star2,prob=0.025),sigsq.ls,D)
avc.pair = ci(AvC,beta.ls,quantile(avc.t.star2,prob=0.975),
	quantile(avc.t.star2,prob=0.025),sigsq.ls,D)
bvc.pair = ci(BvC,beta.ls,quantile(bvc.t.star2,prob=0.975),
	quantile(bvc.t.star2,prob=0.025),sigsq.ls,D)

ci.pair = cbind(avb.pair,avc.pair,bvc.pair)
colnames(ci.pair) = c("A vs. B","A vs. C","B vs. C")
rownames(ci.pair) = c("Lower Bound","Upper Bound")


### latex code
xtable(ci.t)
xtable(ci.res)
xtable(ci.pair)



