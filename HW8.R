####################
## Problem 2
####################
library(graphics)

###MSFT
setwd("C:\\Users\\emurphy\\Dropbox\\MATH352\\HW8")

data = read.table("cholost.txt",header=F)
colnames(data) = c("z","y")
attach(data)

### fit loess
fit1 = loess(y~z,span=0.3,degree=1)
#fit2 = loess(y~z,span=0.3,degree=2)

### plot data and fit
znew = seq(1,100,by=1)
ynew = predict(fit1,newdata=znew)

plot(z,y,cex=0.4)
lines(znew,ynew,col="blue")
#lines(,predict(fit2,newdata=x),col="red")

### bootstrap

B = 50
zboot = matrix(nrow=B,ncol=3)

for(i in 1:B){
	data.star = data[sample(1:nrow(data),replace=T),]
	fit.star = loess(data.star[,2]~data.star[,1],span=0.3,degree=1)
	zboot[i,] = predict(fit.star,c(60,80,100))
}

se = apply(zboot,2,sd)
center = predict(fit1,c(60,80,100))

plot(z,y,cex=0.4)
lines(znew,ynew,col="blue")

arrows(60,center[1]-se[1],60,center[1]+se[1],angle=90,code=3,length=0.1)
arrows(60,center[1]-2*se[1],60,center[1]+2*se[1],angle=90,code=3,length=0.1,col="red")

arrows(80,center[2]-se[2],80,center[2]+se[2],angle=90,code=3,length=0.1)
arrows(80,center[2]-2*se[2],80,center[2]+2*se[2],angle=90,code=3,length=0.1,col="red")

arrows(100,center[3]-se[3],100,center[3]+se[3],angle=90,code=3,length=0.1)
arrows(100,center[3]-2*se[3],100,center[3]+2*se[3],angle=90,code=3,length=0.1,col="red")



