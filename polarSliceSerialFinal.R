###############################################################################
# Serial code for sampling from exp{-|x|}, where |x| is Euclidean norm of vector,
# using polar slice sampler of Roberts/Rosenthal, combining approach presented
# Tibbits Stat Comput (2011) 21:415-430.  Approach uses stepping out procedure
# of Neal:
#	1.  Sample u~U[0,f1x], where f1x = |x|^(d-1) exp(-|x|)
#	2.  Create hypercube s.t. f1(vertices) >= u
#	3.  Sample from density f0x = |x|^-(d-1) over region created in 2.
# Created:  11/06/11 by Edward Murphy 
# Edited: 02/08/12 by Edward Murphy 
#	- Removed exploratory and test code
###############################################################################

install.packages("scatterplot3d")
install.packages("mvtnorm")

library(scatterplot3d)
library(mvtnorm)

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

## initialize
w <- 1      #stepping out parameter (fixed multiple of runif to step out)

x1 <- -10
x2 <- 0
x3 <- 10
X <- rbind ( x1, x2, x3 )

## set parameters of t
df<-1.5				
mu<-0 			## expected value
sigma<-diag(1,nrow=3)  	## covariance matrix

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

			if (max.i > M) M <- max.i

		}

		## sample from mvt until meets AR (u <= f/(M*g) && slice sampler req
		## 	(f1x(cand) > h)

		### NEED TO CHECK IF CANDIDATE IS WITHIN HYPERCUBE??? ##############
		iter<-1
		repeat{
			u<-runif(1)
			Y<-rmvt(n=1,sigma=sigma,df=df)
			if ( u < f0x(t(Y))/(M*dmvt(Y,sigma=sigma,df=df,log=F)) &&
			f1x(t(Y)) > h )		break
			iter<-iter+1
		}

		X <- cbind (X, t(Y)) 
		total.iter<-total.iter+1
	}

return (X[,-1])
}

result<-Sample(X,5000)
par(mfrow=c(3,1))
hist(result[1,],breaks=100)
hist(result[2,],breaks=100)
hist(result[3,],breaks=100)







