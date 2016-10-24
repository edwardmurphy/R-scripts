### Problem 5.4.12
low<--sqrt(5)
up<-sqrt(5)

fx <- function(x) 3/(20*sqrt(5))*(5-x^2)
ex <- function(x) x*3/(20*sqrt(5))*(5-x^2)
ex2 <- function(x) x^2*3/(20*sqrt(5))*(5-x^2)
sig2 <- integrate(ex2,low,up)$value - (integrate(ex,low,up)$value)^2

are.fx <- function(x) (3/(20*sqrt(5))*(5-x^2))^2

are.integrand.val<-integrate(are.fx,low,up)$value

ARE <- 12*sig2*(are.integrand.val)^2
[1] 0.864

### Problem 6.1.4

#data
y<-c(525,570,190,395,370,210,490,250,360,285,630,385,195,295)
x<-c(270,150,270,420,202,255,165,220,305,210,240,300,300,70)

#statistic
z<-y-x

n<-length(z)
k<-n*(n+1)/2

#walsh averages
zstar<-vector("integer",k)
K<-0
for (i in 1:n){
	for (j in i:n){
		K<-K+1
		zstar[K]<-(z[[i]]+z[[j]])/2
	}
}

## functions for cdf of wilcoxen signed rank test
## from Tiffany Head, CGU MATH352
cSRT<-function(n,k){ 
	if(k<0 || n<1)
    	y=0
	else if(n==1)
        if(k==0 || k==1)
        y=1
        else
        y=0
    	else
      y=cSRT(n-1,k-n)+cSRT(n-1,k)
return(y)
}

wilcoxencdf<-function(n,k){
	F=0;
	for (i in 0:k){
    		F=F+cSRT(n,i)
	}
	F=F/(2^n)
	return(F)
}

#d1
1-wilcoxencdf(n,21)
[1] 0.9752808

#CI
upper.ci.limit <- sort(zstar)[k-21]
[1] 232.5

# agree with wilcox.test?
wilcox.test(z,alternative="less",exact=T,correct=F,conf.int=T,
	conf.level=0.975)$conf.int

### Problem 6.1.5

#data
y<-c(11.5,12,9,11.5,13.25,13) #no-exercise
x<-c(9,9.5,9.75,10,13,9.5) #active-exercise

m<-length(y)
n<-length(x)

#statistic
d<-vector("numeric",m*n)
K<-0
for (j in 1:m){
	for (i in 1:n){
		K<-K+1
		d[K] <- y[j]-x[i]
	}
}

## functions for cdf of mann-whitney rank test
## from Tiffany Head, CGU MATH352
tRT<-function(m,n,w){
	if(w<(n*(n+1)/2)) 
    	y=0
	else if(n==0 && w==0)
      y=1
    	else if(m==0 && w==(n*(n+1)/2))
      y=1
      else if(n==0)
      y=0
      else if(m==0)
      y=0
      else
      y=tRT(m,n-1,w-m-n)+tRT(m-1,n,w)
	return(y)
}

manwhitneycdf<-function(m,n,k){
# this form is for the cdf of U = W - n(n+1)/2.  For the cdf of W, set mm=0
# below.
	F=0
	mm=n*(n+1)/2
	for (i in 0:k){
    		F=F+tRT(m,n,i+mm)
	}
	F=F/choose((n+m),n)
	return(F)
}

#d2
manwhitneycdf(6,6,30)
[1] 0.9794372

lower.ci.limit <- sort(d)[m*n+1-30]
[1] -0.5

### Problem 6.2.5

# run data from problem 6.1.4

q = round(n*(n+1)/4 - qnorm(0.975)*sqrt(n*(n+1)*(2*n+1)/24))


