### Compare proposed polar slice sampler using stepping out with
### multivariate t candidate with sampler in Roberts/Rosenthal

######### stepping out/multi t 

library(scatterplot3d)
library(mvtnorm)
library(moments)

################################################################
### f1x function, computes upper limit of vertical slice
### INPUT:  column vector of current MC values
### OUPUT:  vector evaluated at f1x (see intro)

f1x <- function (colvectorX){
	norm.x <- sqrt( t(colvectorX) %*% colvectorX )
	d <- nrow ( as.matrix(colvectorX,ncol=1) )
	eval <- {( norm.x )^( d - 1 )}*{exp ( -norm.x )}
	return ( eval )
}

###
################################################################

################################################################
### f0x function: |x|^{-(d-1)}
### INPUT:  column vector of current MC values
### OUPUT:  vector evaluated at f0x (see intro)

f0x <- function (colvectorX){
	norm.x <- sqrt( t(colvectorX) %*% colvectorX )
	d <- nrow ( as.matrix(colvectorX,ncol=1) )
	eval <- ( norm.x )^(-( d - 1 ))
	return ( eval )
}

###
################################################################

#############################################################################
### matrixVertices function, creates matrix containing all of the vertices
### given the lower and upper bound in each dimension
### INPUTS:  L = vector of lower bound for each dimension
###		 U = vector of upper bound for each dimension
### OUTPUTS:  matrix containing vertices in each column, where nrow 
###		  is dimension size

matrixVertices <- function ( L,U ) {

	## create matrix with lower and upper bounds for each dimension
	LU <- cbind ( L,U ) 

	## start at dimension 1
	A <- matrix(rep(LU[1,],2),nrow=1)

	## loop through remaining dimensions
	for ( i in 2:(nrow(LU)) ) {
		A <- rbind( A, rep(LU[i,],1,each=0.5*ncol(A)))
		if (i != nrow(LU)) A <- cbind(A,A)
	}

	return (A)
} 

###
#############################################################################

#############################################################################
### cubeConstruct function, to create hypercube via stepping out procedure
### takes vector of current MC values and omega (stepping out 
### parameter, and for each dimension first evaluates if f1x at each vertex 
### is < height which is randomly drawn within function; if not, 
### steps out in all dimensions until condition [f(vertex) < height]is met
### INPUTS: height = uniform random variate from step 1 of algorithm,
###		vectorX = vector of current MC values (in rows), 
###		omega = stepping out parameter
### OUTPUTS:  vertices of hypercube s.t. f1x >= height at all vertices

cubeConstruct <- function (height, vectorX, omega){
	u <- runif ( 1 )
	L <- vectorX - omega*u
	U <- L + omega

	vertices <- matrixVertices ( L,U )	
	## easier:  combination function in R??
	
	## evaluate f1x at each vertex
	f.vert <- rep(0,times=ncol(vertices))
	for ( i in 1:ncol(vertices) ){
		f.vert[i] <- f1x(vertices[,i])
	}

	count <- 0 # number of step outs required
	while ( sum(f.vert > height) != 0 ) {  
		count <- count + 1
		L <- L - omega
		U <- U + omega
		vertices <- matrixVertices ( L,U )
		for ( i in 1:ncol(vertices) ){
		f.vert[i] <- f1x(vertices[,i])
		}
	}

	#cat ( "Number of step outs needed:",count,"\n" )
	return ( list(vertices,L,U) )
}

###
#############################################################################

#############################################################################
### accept reject function and optimization to find arg max of f(x)/g(x), where
###	f(x) = |x|^-(p-1), x is p dimensional vector
###	g(x) ~ p-variate t .
### currently need to manually adjust the inputs, write wrapper?
### INPUTS: df = degrees of freedom for mv t
###			changing this flattens the surface (spreads out over support)
###		mu = expected value of mv t
###		sigma = covariance matrix of mv t 
###		vectorX = vector of current MC values (in rows), 
###		omega = stepping out parameter
### OUTPUTS:  vertices of hypercube s.t. f1x >= height at all vertices

ar <- function (x){
	p<-length(x)
	{ t(x) %*% x }^ -( (p-1)/2 )*
	( 1 + 1/df * ( t(x-mu)%*%solve(sigma)%*%(x-mu) ) )^( (df+p)/2 ) 
}

#############################################################################

######  MAIN  #########

Sample<-function(X,n) {
	
	total.iter=1
	while (total.iter <= n) {

		h <- runif ( 1, max = f1x(X[,total.iter]) )
		sample.space <- cubeConstruct ( h, X[,total.iter], w )

		## get M value (for AR; f/g<=M)
		M <- 0
		for (i in 1:ncol(sample.space[[1]])){

			max.i<-(optim(sample.space[[1]][,i],ar, method= "L-BFGS-B", 
			control=list(fnscale=-1),lower=sample.space[[2]]
			,upper=sample.space[[3]]))$value

			if (max.i > M && max.i < 10^5) M <- max.i

		}

		## sample from mvt until meets AR (u <= f/(M*g) && slice sampler req
		## 	(f1x(cand) > h)

		#iter<-1
		repeat{
			u<-runif(1)
			Y<-rmvt(n=1,sigma=sigma,df=df)
			if ( u < f0x(t(Y))/(M*dmvt(Y,sigma=sigma,df=df,log=F)) &&
			f1x(t(Y)) > h &&
			Y[1] > (sample.space[[2]])[[1]] &&
			Y[2] > (sample.space[[2]])[[2]] &&
			Y[3] > (sample.space[[2]])[[3]] &&
			Y[1] < (sample.space[[3]])[[1]] &&
			Y[2] < (sample.space[[3]])[[2]] &&
			Y[3] < (sample.space[[3]])[[3]] )		break
			#iter<-iter+1
		}

		X <- cbind (X, t(Y)) 
		total.iter<-total.iter+1
	}

return (X[,-1])
}

######## Roberts/Rosenthal

PolarSampler <- function(X,n,delta) {

	total.iter=1
	while (total.iter <= n) {
		
		repeat {

			ynew <- runif ( 1, max = f1x(X[,total.iter]) )
			rstar <- log(delta*ynew)/(delta-1)
			rnew <- runif ( 1, max = rstar )
			stdnorm <- rnorm( 3, mean=0, sd=1 )
			theta <- stdnorm / sqrt( sum ( stdnorm^2 ) )
			xnew <- rnew * theta
			if (f1x(xnew) > ynew) {break}
		}

	X<-cbind(X,xnew)
	total.iter<-total.iter+1

	}

	return(X[,-1])
}

## initialize stepping out/multi t

w <- 5      #stepping out parameter (fixed multiple of runif to step out)

#x1 <- -10
#x2 <- 0
#x3 <- 10
x1 <- runif(1,-10,10)
x2 <- runif(1,-10,10)
x3 <- runif(1,-10,10)
Xstart <- rbind ( x1, x2, x3 )

## set parameters of t
df<-1.5				
mu<-0 			## expected value
sigma<-diag(1,nrow=3)  	## covariance matrix
#sigma<-diag(10,nrow=3)
#sigma<-matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),nrow=3,ncol=3)
#sigma<-matrix(c(10,1,1,1,10,1,1,1,10),nrow=3,ncol=3)
#sigma<-diag(0.7,nrow=3) 

## initialize Roberts/Rosenthal

delta <- 0.05
#delta <- 0.0375

## run both

system.time(result<-Sample(Xstart,10000))
system.time(polar<-PolarSampler(Xstart,10000,delta))

### compare moments
apply(result,1,mean)
apply(polar,1,mean)
apply(result,1,var)
apply(polar,1,var)
apply(result,1,kurtosis)
apply(polar,1,kurtosis)

min(apply(result,1,var))
max(apply(result,1,var))
min(apply(polar,1,var))
max(apply(polar,1,var))

min(apply(result,1,kurtosis))
max(apply(result,1,kurtosis))
min(apply(polar,1,kurtosis))
max(apply(polar,1,kurtosis))

### graph results
par(mfrow=c(2,3))

hist(result[1,],breaks=50,xlim=c(-10,10))
hist(result[2,],breaks=50,xlim=c(-10,10))
hist(result[3,],breaks=50,xlim=c(-10,10))
hist(polar[1,],breaks=50,xlim=c(-10,10))
hist(polar[2,],breaks=50,xlim=c(-10,10))
hist(polar[3,],breaks=50,xlim=c(-10,10))

### effect of increasing stepping out parameter w
### fixed covariance structure

windows(record=TRUE)
pdf("stepping out hists sigma I.pdf")
w <- 1
result<-Sample(Xstart,10000)
sampler.var<-apply(result,1,var)

par(mfrow=c(2,3))
hist(result[1,],breaks=50,xlim=c(-10,10),main=paste("x1, w=",w))
hist(result[2,],breaks=50,xlim=c(-10,10),main=paste("x2, w=",w))
hist(result[3,],breaks=50,xlim=c(-10,10),main=paste("x3, w=",w))
hist(polar[1,],breaks=50,xlim=c(-10,10),main=paste("x1, Roberts/Rosenthal"))
hist(polar[2,],breaks=50,xlim=c(-10,10),main=paste("x1, Roberts/Rosenthal"))
hist(polar[3,],breaks=50,xlim=c(-10,10),main=paste("x1, Roberts/Rosenthal"))

repeat {
	w<-w+1
	result<-Sample(Xstart,10000)
	sampler.var<-cbind(sampler.var,apply(result,1,var))

	hist(result[1,],breaks=50,xlim=c(-10,10),main=paste("x1, w=",w))
	hist(result[2,],breaks=50,xlim=c(-10,10),main=paste("x2, w=",w))
	hist(result[3,],breaks=50,xlim=c(-10,10),main=paste("x3, w=",w))
	hist(polar[1,],breaks=50,xlim=c(-10,10),main=paste("x1, Roberts/Rosenthal"))
	hist(polar[2,],breaks=50,xlim=c(-10,10),main=paste("x1, Roberts/Rosenthal"))
	hist(polar[3,],breaks=50,xlim=c(-10,10),main=paste("x1, Roberts/Rosenthal"))

	if (w == 10) break
}

dev.off()

pdf("var plot sigma I.pdf")
plot(c(1:10),sampler.var[1,],xlab="Stepping Out Parameter",ylab="Variance",
	main="Variance by Stepping Out Parameter, Sigma = Identity")
dev.off()


	
	
	



