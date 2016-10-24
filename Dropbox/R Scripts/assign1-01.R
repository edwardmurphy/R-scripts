alph = c(1,2,3,4)
lam = c(1,10,100,200)

x = matrix(0,nrow=length(alph),ncol=length(lam))
for (i in 1:length(alph)){
	for (j in 1:length(lam)){
		x[i,j] = alph[i]*lam[j]
	}
}


true = matrix(-1,nrow=length(alph),ncol=length(lam))
for (i in 1:length(alph)){
	for (j in 1:length(lam)){
		true[i,j] = dpois(x[i,j],lam[j])
	}
}

approx <- matrix(-1,nrow=length(alph),ncol=length(lam))
for (i in 1:length(alph)){
	for (j in 1:length(lam)){
		approx[i,j] = 1/sqrt(2*pi*alph[i]*lam[j])*(exp(alph[i])/(exp(1)*alph[i]^alph[i]))^lam[j]
	}
}

approx_alph1 <- matrix(-1,nrow=1,ncol=length(lam))
for (j in 1:length(lam)){
		approx_alph1[1,j] = 1/sqrt(2*pi*lam[j])
}

true
approx
approx_alph1


library(xtable)
xtable(true,display=c("d","E","E","E","E"))

xtable(approx,display=c("d","E","E","E","E"))

xtable(approx_alph1,display=c("d","E","E","E","E"))
