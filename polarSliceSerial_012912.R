install.packages("scatterplot3d")
library(scatterplot3d)

f1x <- function (vectorX){
	norm.x <- norm ( as.matrix(vectorX,ncol=1), type = "f" )
	d <- nrow ( as.matrix(vectorX,ncol=1) )
	eval <- {( norm.x )^( d - 1 )}*{exp ( -norm.x )}
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

## initialize
w <- 1      #stepping out parameter (fixed multiple of runif to step out)


x1 <- -5
x2 <- 0
x3 <- 5
X <- rbind ( x1, x2, x3 )

h <- runif ( 1, max = f1x(X) )
sample.space <- cubeConstruct ( h, X, w )


df<-5				
mu<-0 			
sigma<-diag(1,nrow=3) 

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

optim(c(-1,0,1),ar, method= "L-BFGS-B", 
	control=list(fnscale=-1),lower=c(sample.space[[2]])
	,upper=sample.space[[3]])



 

