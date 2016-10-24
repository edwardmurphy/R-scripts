setwd("C:\\Users\\Owner\\Dropbox\\MATH353")
norm_rand_MC <- rnorm(1000000)
x <- seq(-5,5,by=0.0001)
y <- dnorm(x)

pdf("hists_assign3.pdf")
par(mfrow=c(2,1))
hist(norm_rand_MC,breaks=50,freq=F,main=
"Histogram of Std Norm Variables via Mersenne-Twister 
with True Density Overlay")
points(x,y,type="l",lty=2)


iter<-1000000
z<-vector(mode="numeric",length=iter)
for(i in 1:iter){
	z[i]<-sum(runif(12))-6
}
	
hist(z,breaks=50,freq=F,main="Histogram of Std Norm Variables
from Uniform (n=12) with True Density Overlay")
points(x,y,type="l",lty=2)
dev.off()


C <- sort(z)
y <- seq(1,iter)/(iter+1)
y_true <- pnorm(x)
png("cdf_assign3.png")
par(mfrow=c(1,1))
plot(C,y,type="l",xlab="x",ylab="F(x)")
points(x,y_true,type="l",lty=2)
legend("topleft",legend=c("Empirical CDF","Normal CDF"),lty=c(1,2))
dev.off()

