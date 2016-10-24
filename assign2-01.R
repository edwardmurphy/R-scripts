setwd("/home/ed/Dropbox/MATH353")
y<-seq(-10,10,by=0.001)

#CDFs
lim_cdf <- exp(-exp(-y))
n1_cdf <- (1-(exp(-y)/1))^1
n2_cdf <- (1-(exp(-y)/2))^2
n5_cdf <- (1-(exp(-y)/5))^5

#PDFs
lim_pdf <- exp(-y)*exp(-exp(-y))
n1_pdf <- 1*(1-(exp(-y)/1))^(1-1)*exp(-y)
n2_pdf <- 2*(1-(exp(-y)/2))^(2-1)*exp(-y)
n5_pdf <- 5*(1-(exp(-y)/5))^(5-1)*exp(-y)

pdf("assign2plots.pdf")
par(mfrow=c(3,1))

#CDFs
plot(y,lim_cdf,lty=1,type="l",main="CDFs",ylab="Probability")
points(y,n1_cdf,lty=2,type="l")
points(y,n2_cdf,lty=3,type="l")
points(y,n5_cdf,lty=4,type="l")
legend("topleft",legend=c("limit","n=1","n=2","n=5"),lty=c(1,2,3,4))

#PDFs limit scale
plot(y,lim_pdf,lty=1,type="l",main="PDFs, limit scale",ylab="Density")
points(y,n1_pdf,lty=2,type="l")
points(y,n2_pdf,lty=3,type="l")
points(y,n5_pdf,lty=4,type="l")
legend("topleft",legend=c("limit","n=1","n=2","n=5"),lty=c(1,2,3,4))

#PDFs divergent scale
plot(y,lim_pdf,lty=1,type="l",ylim=c(0,10),main="PDFs, divergent scale",ylab="Density")
points(y,n1_pdf,lty=2,type="l")
points(y,n2_pdf,lty=3,type="l")
points(y,n5_pdf,lty=4,type="l")
legend("topleft",legend=c("limit","n=1","n=2","n=5"),lty=c(1,2,3,4))

dev.off()


