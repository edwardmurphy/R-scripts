### M/L polar slice sampler, generalized for any target 
### and any dimension

library(mvtnorm)
library(moments)

################################################################
### target density function
### INPUT:  column vector of current MC values
### OUPUT:  vector evaluated at target (see intro)

target <- function (colvectorX){

	norm.x <- sqrt( t(colvectorX) %*% colvectorX )
	eval <- (exp ( -norm.x ))
	return ( eval )
}

###
################################################################


################################################################
### f1x function, computes upper limit of vertical slice
###	equal to |x|^(d - 1) * target
###	so f0x * f1x = target
### INPUT:  column vector of current MC values
### OUPUT:  vector evaluated at f1x (see intro)

f1x <- function (colvectorX){

	norm.x <- sqrt( t(colvectorX) %*% colvectorX )
	d <- nrow ( as.matrix(colvectorX,ncol=1) )
	eval <- {( norm.x )^( d - 1 )}*target(colvectorX)
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

	#count <- 0 # number of step outs required
	while ( sum(f.vert > height) != 0 ) {  
		#count <- count + 1
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
	f1x.obs<-f1x(X[,total.iter])
	h.obs<-0
	M.obs<-0
	iter.obs<-0
	
	height.time <- 0
	sample.space.time <- 0
	M.time <- 0
	AR.time <- 0

	repeat {

		height.time[total.iter]<- system.time(
			h <- runif ( 1, max = f1x(X[,total.iter]) ) )[3]
		sample.space.time[total.iter] <- system.time(
			sample.space <- cubeConstruct ( h, X[,total.iter], w ) )[3]

		## get M value (for AR; f/g<=M)
		M <- 0
		M.time[total.iter]<-system.time(
			for (i in 1:ncol(sample.space[[1]])){

				max.i<-(optim(sample.space[[1]][,i],ar, method= "L-BFGS-B", 
				control=list(fnscale=-1),lower=sample.space[[2]]
				,upper=sample.space[[3]]))$value

				if (max.i > M && max.i < 10^5) M <- max.i

			} )[3]

		## sample from mvt until meets AR (u <= f/(M*g) && slice sampler req
		## 	(f1x(cand) > h)

		iter<-1
		AR.time[total.iter]<-system.time(
			repeat{
				u<-runif(1)
				Y<-rmvt(n=1,sigma=sigma,df=df)

				lower<-rep(0,length(Y))
				upper<-rep(0,length(Y))
				for (i in 1:length(Y)){
					if (Y[i] > sample.space[[2]][[i]]) lower[i]<-1
					if (Y[i] < sample.space[[3]][[i]]) upper[i]<-1
				}

				if ( u < f0x(t(Y))/(M*dmvt(Y,sigma=sigma,df=df,log=F)) &&
				f1x(t(Y)) > h &&
				sum(lower) == length(Y) &&
				sum(upper) == length(Y)) break
				iter<-iter+1
			} )[3]
		
		X <- cbind (X, t(Y)) 
		f1x.obs<-rbind(f1x.obs,f1x(t(Y)))
		h.obs<-rbind(h.obs,h)
		M.obs<-rbind(M.obs,M)
		iter.obs<-rbind(iter.obs,iter)
		
		total.iter<-total.iter+1
		if (total.iter == (n + 1)) break
	}

	# return all accepted X values (minus initialization), the last h value (for
	#	investigation if break), the last sample.space (investigation), 
	#	the last M (investigation), the last u (investigation), the last 
	#	Y (investigation), all f1x densities at accepted values (minus 
	#	initialization), all random heights that lead to accepted values (minus
	#	initialization), all M values that lead to accepted values (minus
	#	initialization), all iter values (number of draws for accept-reject) that 
	#	lead to accepted values (minus initialization), and timed steps (height, 
	#	sample.space, M optimization, AR for each iteration)

	return (list(X[,-1], h, sample.space, M, u, Y, f1x.obs[-1,], 
		h.obs[-1,], M.obs[-1,], iter.obs[-1,], height.time, sample.space.time,
		M.time, AR.time))
}

######## Roberts/Rosenthal

PolarSampler <- function(X,n,delta) {

	total.iter=1
	while (total.iter <= n) {
		
		repeat {

			ynew <- runif ( 1, max = f1x(X[,total.iter]) )
			rstar <- log(delta*ynew)/(delta-1)
			rnew <- runif ( 1, max = rstar )
			stdnorm <- rnorm( 10, mean=0, sd=1 )
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

w <- 1     #stepping out parameter (fixed multiple of runif to step out)

Xstart <- t(t(runif(10,-10,10)))

## set parameters of t
df<-1.5				
mu<-0 
			## expected value
sigma<-diag(1,nrow=length(Xstart))  	## covariance matrix

## initialize Roberts/Rosenthal

#delta <- 0.05
delta <- 0.0375

## run both

time<-system.time(result<-Sample(Xstart,20)) 
system.time(polar<-PolarSampler(Xstart,250,delta))

# for dim=10, 2500 samples takes 8.3 hours with M/L


## get results from list 
X <- result[[1]] 
h <- result[[2]]
sample.space <- result[[3]]
M <- result[[4]]
u <- result[[5]]
Y <- result[[6]]
f1x.obs <- result[[7]]
h.obs <- result[[8]]
M.obs <- result[[9]]
iter.obs <- result[[10]]
height.time <- result[[11]]
sample.space.time <- result[[12]]
M.time <- result[[13]]			#90% of the total!!
AR.time <- result[[14]]

iter.time<-height.time+sample.space.time+M.time+AR.time


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



	
	
	



