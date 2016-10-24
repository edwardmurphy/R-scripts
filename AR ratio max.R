###### some AR maximization fun

## get a hypercube

target <- function (colvectorX){

	norm.x <- sqrt( t(colvectorX) %*% colvectorX )
	eval <- (exp ( -norm.x ))
	return ( eval )
}

f1x <- function (colvectorX){

	norm.x <- sqrt( t(colvectorX) %*% colvectorX )
	d <- nrow ( as.matrix(colvectorX,ncol=1) )
	eval <- {( norm.x )^( d - 1 )}*target(colvectorX)
	return ( eval )
}

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

ar <- function (x){
	p<-length(x)
	{ t(x) %*% x }^ -( (p-1)/2 )*
	( 1 + 1/df * ( t(x-mu)%*%solve(sigma)%*%(x-mu) ) )^( (df+p)/2 ) 
}

X<-runif(10,-10,10)
w<-5
df<-1.5
mu<-0
sigma<-diag(1,10)

h <- runif ( 1, max = f1x(X) )
sample.space <- cubeConstruct ( h, X, w ) 

# get list of M values starting from each vertex
vertex.M<-0
for (i in 1:ncol(sample.space[[1]])){

	vertex.M[i]<-optim(sample.space[[1]][,i],ar, method= "L-BFGS-B", 
		control=list(fnscale=-1),lower=sample.space[[2]]
		,upper=sample.space[[3]])$value
}

plot(vertex.M,pch=1)

# get list of ar ratio evaluated at each vertex

vertex.dens<-0
for (i in 1:ncol(sample.space[[1]])){

	vertex.dens[i]<-ar(sample.space[[1]][,i])
}

points(vertex.dens,col="red",pch=2)

## optim function, change initial values- is vertex the arg max?

optim(sample.space[[1]][,90],ar, method= "L-BFGS-B", 
		control=list(fnscale=-1),lower=sample.space[[2]]
		,upper=sample.space[[3]])

## code a loop to randomly draw initial values, find max
## finds the vertex from random initialization

par<-rep(0,10)
m<-0

for (i in 1:1000){
	x<-runif(10,sample.space[[2]],sample.space[[3]])
	max<-optim(x,ar, method= "L-BFGS-B", 
		control=list(fnscale=-1),lower=sample.space[[2]]
		,upper=sample.space[[3]])$value
		
	if (max > m) {
		m<-max
		par<-optim(x,ar, method= "L-BFGS-B", 
		control=list(fnscale=-1),lower=sample.space[[2]]
		,upper=sample.space[[3]])$par
	}
}

}


x1<-0
x2<-0
x3<-0
x4<-0
x5<-0
y<-1
for (i in 1:10000){
	x<-runif(10,-10,10)
	x1[i]<-x[1]
	x2[i]<-x[2]
	x3[i]<-x[3]
	x4[i]<-x[4]
	x5[i]<-x[5]
	y[i]<-ar(x)
}

scatterplot3d(x1,x2,y,box=F, angle=24, pch=20)





