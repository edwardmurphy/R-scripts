RSD <- 0
sigma <- 0
RSDcorrected <- 0  # correction for small to med sample sizes
sample.for.rsd <- function (m,n,mu,rsd){
	
	for (i in 1:m){
		iter <- rnorm(n,mean=mu,sd=mu*rsd/100)
		sigma[i] <-  sd(iter)	
		RSD[i] <- sd(iter)/mean(iter)*100
		RSDcorrected[i] <- (1 + 1/(4*n))*RSD[i]
	}
	return(list(sigma,RSD,RSDcorrected))
}

n3 <- sample.for.rsd(50000,3,106,2)
n5 <- sample.for.rsd(50000,5,106,2)
n50 <- sample.for.rsd(50000,50,106,2)

# how std deviation changes with sample size
summary(n3[[1]])
summary(n5[[1]])
summary(n50[[1]])

# how RSD changes with sample size
summary(n3[[2]])
summary(n5[[2]])
summary(n50[[2]])

# how RSD corrected changes with sample size
summary(n3[[3]])
summary(n5[[3]])
summary(n50[[3]])

par(mfrow=c(3,1))
hist(n3[[1]],breaks=100)
hist(n5[[1]],breaks=100)
hist(n50[[1]],breaks=100)

a <- 2
for (i in 1:length(n3[[1]])){
	if (n3[[1]][i]<=0.046) a[i] <- 1
	else a[i] <- 0
}

sum(a)/length(a)






