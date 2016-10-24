### R code from vignette source 'HW3.Rnw'

###################################################
### code chunk number 1: HW3.Rnw:16-49
###################################################
##### parameters
lam = 0.001
n = 4
k = 2
tau = 500

##### MTTF without renewal
sum <- 0
for (i in k:n){
	sum <- sum + 1/i
}

mttf <- 1/lam*sum

##### MTTF with renewal
num<-0
for(i in k:n){
	for(j in 0:(n-i)){
		num <- num + {gamma(n+1)/{gamma(i+1)*gamma(j+1)*gamma(n-i-j+1)}* 			(-1)^(j+1) * 1/{lam*(i+j)} * {exp(-((i+j)*lam*tau)) - 1}}
	}
}

den<-0
for(i in k:n){
	for(j in 0:(n-i)){
		den <- den + {gamma(n+1)/{gamma(i+1)*gamma(j+1)*gamma(n-i-j+1)}* 			(-1)^j * exp(-((i+j)*lam*tau))}
	}
}

den <- 1-den
mttfPM <- num/den

inc<-mttfPM/mttf


###################################################
### code chunk number 2: rel_compare
###################################################
##### reliability function without renewal
R_t <- function(lam,n,k,t){
	f<-0
	for(i in k:n){
		for(j in 0:(n-i)){
			f <- f + {gamma(n+1)/{gamma(i+1)*gamma(j+1)*gamma(n-i-j+1)}* 				(-1)^j * exp(-((i+j)*lam*t))}
		}
	}
	return(f)
}

##### reliability function with renewal
R_tau <- function(lam,n,k,t,tau){
	f<-{R_t(lam,n,k,tau)}^(floor(t/tau))*R_t(lam,n,k,t-tau*floor(t/tau))
	return(f)
}

##### plot reliabilities
t <- seq(0,4000,by=0.5)
f <- R_t(lam,n,k,t)
g <- R_tau(lam,n,k,t,tau)

png("rel_compare.png")
plot(t,f,type="l",xlab="Time",ylab="Reliability")
par(new=TRUE)
points(t,g,type="l",col="blue")
legend("topright",c("No PM","PM"),lty=c(1,1),col=c("black","blue"))
null<-dev.off()


