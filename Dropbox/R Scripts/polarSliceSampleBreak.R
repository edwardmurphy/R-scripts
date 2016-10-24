### Performs polar slice with stepping out using a break in sampler code
### to investigate the cause of sampler "hold-up" in most cases
### stores and prints critical results from sampler

### 03/07/12:  cause of hold-up appears to be large M values
### 			perhaps when near 0 for domain????

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

Sample<-function(X,n) {
	
	total.iter=1
	f1x.obs<-f1x(X[,total.iter])
	h.obs<-0
	M.obs<-0
	iter.obs<-0
	repeat {

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

		iter<-1
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
			iter<-iter+1
			if (iter > 29999) cat("Iteration has reached 30000 \n")
			if (iter > 29999) break
		}
		# first tried 10,000 but this caused sampler to break when h was near the 
		# max of f1x, so in a location where number of samples expected to be high
		# but normal part of algorithm
		if (iter > 29999) break

		## only perform the following if accept candidate
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
	#	lead to accepted values (minus initialization)

	return (list(X[,-1], h, sample.space, M, u, Y, f1x.obs[-1,], 
		h.obs[-1,], M.obs[-1,], iter.obs[-1,]))
}

## initialize stepping out/multi t

w <- 1      #stepping out parameter (fixed multiple of runif to step out)

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

#### RUN SAMPLER

system.time(result<-Sample(Xstart,10000))

####

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

# find vertex where M is large

for (i in 1:ncol(sample.space[[1]])){

	cat("Vertex",i,", M =",(optim(sample.space[[1]][,i],ar, method= "L-BFGS-B", 
		control=list(fnscale=-1),lower=sample.space[[2]]
		,upper=sample.space[[3]])$value),"\n")
}

## replace i in  sample.space[[1]][,i]  in call below with vertex from for loop above
## to view the arg max. 

optim(sample.space[[1]][,1],ar, method= "L-BFGS-B", 
		control=list(fnscale=-1),lower=sample.space[[2]]
		,upper=sample.space[[3]])

max(M.obs)	# for sampler that reached 10,000 iterations, this was 4568



## find max of f1x
## coded when thought cause of hold-up was h very near f1x max

f1x.optim <- function (x){
	p <- length(x)
	{( sqrt( t(x) %*% x ) )^( p - 1 )}*{exp ( -(sqrt( t(x) %*% x )) )}
}

optim(c(-10,0,10),f1x.optim, method= "Nelder-Mead", 
			control=list(fnscale=-1))


