###############################################################################
# Serial code for sampling from exp{-|x|}, where |x| is Euclidean norm of vector,
# using polar slice sampler of Roberts/Rosenthal, combining approach presented
# Tibbits Stat Comput (2011) 21:415-430.  Approach uses stepping out procedure
# of Neal:
#	1.  Sample u~U[0,f1x], where f1x = |x|^(d-1) exp(-|x|)
#	2.  Create hypercube s.t. f1(vertices) >= u
#	3.  Sample from density f0x = |x|^-(d-1) over region created in 2.
# Created:  11/06/11 by Edward Murphy 
###############################################################################

install.packages("scatterplot3d")
install.packages("mvtnorm")

library(scatterplot3d)
library(mvtnorm)

################################################################
### f1x function, computes upper limit of vertical slice
### INPUT:  vector of current MC values
### OUPUT:  vector evaluated at f1x (see intro)

f1x <- function (vectorX){
	norm.x <- sqrt( t(vectorX) %*% vectorX )
	d <- nrow ( as.matrix(vectorX,ncol=1) )
	eval <- {( norm.x )^( d - 1 )}*{exp ( -norm.x )}
	return ( eval )
}

# test 
x1 <- 1.2
x2 <- 2.3
x3 <- 3.6
X <- rbind ( x1, x2, x3 )
f1x(X)  #by calculator, 0.23288073087

###
################################################################

################################################################
### f0x function: norm(x)^-(d-1)
### INPUT:  vector of current MC values
### OUPUT:  vector evaluated at f0x (see intro)

f0x <- function (vectorX){
	norm.x <- sqrt( t(vectorX) %*% vectorX )
	d <- nrow ( as.matrix(vectorX,ncol=1) )
	eval <- ( norm.x )^(-( d - 1 ))
	return ( eval )
}

# test 
x1 <- 1.2
x2 <- 2.3
x3 <- 3.6
X <- rbind ( x1, x2, x3 )
f0x(X)  #by calculator, 0.0507872016

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
	while ( sum(f.vert > h) != 0 ) {  
		count <- count + 1
		L <- L - omega
		U <- U + omega
		vertices <- matrixVertices ( L,U )
		for ( i in 1:ncol(vertices) ){
		f.vert[i] <- f1x(vertices[,i])
		}
	}

	cat ( "Number of step outs needed:",count,"\n" )
	return ( list(vertices,L,U) )
}

###
#############################################################################

######  MAIN  #########

## initialize
w <- 1      #stepping out parameter (fixed multiple of runif to step out)
# w <- 3

x1 <- -1
# x1 <- -0.5
x2 <- 0
# x2 <- 0
x3 <- 1
# x3 <- 0.5
X <- rbind ( x1, x2, x3 )
## does not work at all 0's, since h is then automatically zero

h <- runif ( 1, max = f1x(X) )
sample.space <- cubeConstruct ( h, X, w )

scatterplot3d((sample.space[[1]])[1,],(sample.space[[1]])[2,],(sample.space[[1]])[3,], 
	box=F, angle=24, pch=20 )

#############################################################################
### accept reject function and optimization to find arg max of f(x)/g(x), where
###	f(x) = |x|^(-p-1), x is p dimensional vector
###	g(x) ~ p-variate t .
### currently need to manually adjust the inputs, write wrapper?
### INPUTS: df = degrees of freedom for mv t
###			changing this flattens the surface (spreads out over support)
###		mu = expected value of mv t
###		sigma = covariance matrix of mv t 
###		vectorX = vector of current MC values (in rows), 
###		omega = stepping out parameter
### OUTPUTS:  vertices of hypercube s.t. f1x >= height at all vertices

## set parameters of t

df<-1.5				
mu<-0 			## expected value
sigma<-diag(1,nrow=3)  	## covariance matrix

####### EXPLORATORY STUFF
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

### test optim function

upper.ellipsoid<-function(x){
	x1<-x[1]
	x2<-x[2]
	sqrt( 16-4*x1^2-x2^2 )
}

optim(c(-1,1),upper.ellipsoid,control=list(fnscale=-1)) #true: x1=0,x2=0,@4)
#######

### accept reject optimization

ar <- function (x){
	p<-length(x)
	x1<-x[1]
	x2<-x[2]
	x3<-x[3]
	{x1^2 + x2^2 + x3^2}^ -( (p-1)/2 )*
	( 1 + 1/df * ( t(x-mu)%*%solve(sigma)%*%(x-mu) ) )^( (df+p)/2 ) 
}
	
optim(c(-1,0,1),ar, control=list(fnscale=-1))  # uh oh again

optim(c(0,-1,1),ar, method= "L-BFGS-B", 
	control=list(fnscale=-1),lower=c(sample.space[[2]])
	,upper=sample.space[[3]])

## likely to be local maxima, so choice of initial values is key
## use the points in sample.space[[1]] (so check 2^d points?)

optim(sample.space[[1]][,1],ar, method= "L-BFGS-B", 
	control=list(fnscale=-1),lower=c(sample.space[[2]])
	,upper=sample.space[[3]])


M <- 0
for (i in 1:ncol(sample.space[[1]])){

	max.i<-(optim(sample.space[[1]][,i],ar, method= "L-BFGS-B", 
	control=list(fnscale=-1),lower=c(sample.space[[2]])
	,upper=sample.space[[3]]))$value

	if (max.i > M) M <- max.i

}

## sample from mvt until meets AR (u <= f/(M*g) & meets slice sampler req
iter<-1
repeat{
	u<-runif(1)
	Y<-rmvt(n=1,sigma=sigma,df=df)
	if ( u < f0x(t(Y))/(M*dmvt(Y,sigma=sigma,df=df,log=F)) &&
	  f1x(t(Y)) > h )		break
	iter<-iter+1
}
iter
Y






###
#############################################################################






#######################################################################
### need to sample from |x|^-(d-1) within sample.space
### following code is exploratory, trying to get a sense of what target 
### density looks like for R^3 vector

### sample randomly within hypercube
x1.rand <- runif(100000,sample.space[1,1],sample.space[1,2])
x2.rand <- runif(100000,sample.space[2,1],sample.space[2,3])
x3.rand <- runif(100000,sample.space[3,1],sample.space[3,5])

scatterplot3d(x1.rand,x2.rand,x3.rand, box=F, angle=24, pch=20)

f0x <- { (x1.rand^2 + x2.rand^2 + x3.rand^2) ^ (1/2) } ^ (-2)
### code this as function

par(mfrow=c(3,1))
plot(x1.rand,f0x, pch=20, ylim=c(0,1))
plot(x2.rand,f0x, pch=20, ylim=c(0,1))
plot(x3.rand,f0x, pch=20, ylim=c(0,1))

par(mfrow=c(1,1))
scatterplot3d(x1.rand,x2.rand,f0x, box=F, angle=24, pch=20, zlim=c(0.001,1))
scatterplot3d(x1.rand,x2.rand,f0x, box=F, angle=24, pch=20, zlim=c(0.01,1))
scatterplot3d(x1.rand,x2.rand,f0x, box=F, angle=24, pch=20, zlim=c(0.1,1))

### nearly uniform unless near/at zero?? 
### Q: Is peak always at zero?? PROVE THIS!!!!!!!!

### now we have cube (value to next sampling step unknown), need to sample from 
### |x|^-(d-1) within this cube
### use accept-reject and maximize to find M
### candidate: multivariate t

### what does actual target density look like?



fx <- exp(-(x1.rand^2+x2.rand^2+x3.rand^2)^0.5)

par(mfrow=c(3,1))
plot(x1.rand,fx)
plot(x2.rand,fx)
plot(x3.rand,fx)

par(mfrow=c(1,1))
scatterplot3d(x1.rand,x2.rand,fx)

### let's sample uniformly from cube, see what happens!

f1x.sample <- h/2
samp <- 0
while (sum(f1x.sample > h) == 0) {
	samp <- samp + 1
	x1.sample <- runif(1,sample.space[1,1],sample.space[1,2])	
	x2.sample <- runif(1,sample.space[2,1],sample.space[2,3])
	x3.sample <- runif(1,sample.space[3,1],sample.space[3,5])

	X.sample <- rbind(x1.sample, x2.sample, x3.sample)

	f1x.sample <- f1x(X.sample)
}

cat("Number of sampling iterations:",samp,"\n")
#######################################################################



















	






	