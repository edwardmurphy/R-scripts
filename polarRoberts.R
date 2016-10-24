f1x <- function (colvectorX){
	norm.x <- sqrt( t(colvectorX) %*% colvectorX )
	d <- nrow ( as.matrix(colvectorX,ncol=1) )
	eval <- {( norm.x )^( d - 1 )}*{exp ( -norm.x )}
	return ( eval )
}



## Initialize

delta <- 0.05

x1 <- -10
x2 <- 0
x3 <- 10
X <- rbind ( x1, x2, x3 )


PolarSampler <- function(X,n,delta) {

	total.iter=1
	while (total.iter <= n) {
		
		repeat {

			ynew <- runif ( 1, max = f1x(X[,total.iter]) )
			rstar <- -log(delta*ynew)/(2*delta)
			rnew <- runif ( 1, max = rstar )
			stdnorm <- rnorm( 3, mean=0, sd=1 )
			theta <- stdnorm / sqrt( sum ( stdnorm^2 ) )
			xnew <- rnew * theta
			if (f1x(xnew) > ynew) {break}
		}

	X<-cbind(X,xnew)
	total.iter<-total.iter+1

	}

	return(X)
}
		