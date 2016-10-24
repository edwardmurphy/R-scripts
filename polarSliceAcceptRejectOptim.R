install.packages("scatterplot3d")
library(scatterplot3d)

#### Code to find arg max of f(x)/g(x), where
####	f(x) = |x|^(-p-1), x is p dimensional vector
####	g(x) ~ p-variate t 

## set parameters of t

df<-5			## degrees of freedom, changing this flattens the
				## parabolic shape (higher df, more compressed)
mu<-0 			## expected value
sigma<-diag(1,nrow=3)  	## covariance matrix

## first, let's graph the ratio

# generate random support

x1<-runif(10000,-5,5)
x2<-runif(10000,-5,5)
x3<-runif(10000,-5,5)
supp.rand<-cbind(x1,x2,x3)

p<-dim(supp.rand)[[2]]

# evaluate at support

fx<-rep(0,length(x1))
for (i in 1:nrow(supp.rand)){
	x<-rbind(x1[i],x2[i],x3[i])   	# 3x1
	MU<-cbind(rep(mu,length(x)))	# 3x1

	fx[i]<- norm(x)^(-(p-1)) *
		( 1 + 1/df * ( t((x-MU)) %*% solve(sigma) %*% (x-MU) ) )^( (df+p)/2 ) 
}

scatterplot3d(x1,x2,fx,box=F, angle=24, pch=20)		#uh oh
scatterplot3d(x1,x3,fx,box=F, angle=24, pch=20)
scatterplot3d(x2,x3,fx,box=F, angle=24, pch=11)

## now use optim function to find arg max and max

# test function first

upper.ellipsoid<-function(x){
	x1<-x[1]
	x2<-x[2]
	sqrt( 16-4*x1^2-x2^2 )
}

optim(c(-1,1),upper.ellipsoid,control=list(fnscale=-1)) #true: x1=0,x2=0,@4)

# now accept-reject ratio

ar <- function (x){
	p<-length(x)
	x1<-x[1]
	x2<-x[2]
	x3<-x[3]
	{x1^2 + x2^2 + x3^2}^ -( (p-1)/2 )*
	( 1 + 1/df * ( t(x-mu)%*%solve(sigma)%*%(x-mu) ) )^( (df+p)/2 ) 
}
	
optim(c(-1,0,1),ar, control=list(fnscale=-1))  # uh oh again



