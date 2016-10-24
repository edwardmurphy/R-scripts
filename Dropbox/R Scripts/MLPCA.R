## initial attempts at MLPCA- NOT OPTIMIZED, VERY SLOW
## use MLPCAfxn for optimized version
## performed some optimization on 03.19.13

## read data 
data<-read.csv("U://data.csv",header=F,sep=",")
data<-read.csv("/home/ed/Dropbox/Bayesian PCA/RawData.csv",header=F,sep=",")
data<-read.csv("C:\\Users\\emurphy\\Dropbox\\Bayesian PCA\\RawData.csv",
	header=F,sep=",")

## maturities
mat<-c(0.5,1,2,3,5,7,10)

#### exploratory

  plot(mat,data[300,],type="l",
  #ylim=c(0,.17)
  )

#### hyperparameters of MLPCA function

  ## set size of truncation (number of components)
  p<-3 #typical of bond problem

  ## start super easy with covariances
  
    psi<-diag(0.001,nrow=nrow(data))
    ## same var for each month for any maturity
    ## value from sample variance of each column in dataset
    ## appear to all be ~0.001
    # psi <- diag(0.1, nrow=nrow(data))
    
    sigma<-diag(0.00006, nrow=ncol(data))
    ## same var for each maturity for any month
    # sigma <- diag(0.1, nrow=ncol(data))
 
  ## convergence rates 
  delta<-1
  conv<-10^-10
  
    
### the following are sample covariance matrices
### not used yet since covariance structure is faked 
### for algorithmic learning reasons
### still to determine role of sample covariance matrix
  
  ## cov matrix of column space
  cov.col<-cov(data)

  ## cov matrix of row space
  cov.row<-cov(t(data))

##### ALTERNATING LEAST SQUARES FOR MLPCA #####  
  


  MLPCA<-function(X,psi,sigma,p,delta,conv){

    ## "column space" := original dimension of X  
    m = nrow(X)
    n = ncol(X)

    ## transpose data
    ## call this the row space (use columns of tranpose)
    Xprime<-t(X)

    ## for now only need to compute inverses of psi,sigma once
    ## will change if different variances later
    psi.inv <- solve(psi)
    sig.inv <- solve(sigma)
   
    ## create matrices for max likelihood estimates in col and row space
    X.maxlike.col<-matrix(0,nrow=m,ncol=n)
    X.maxlike.row<-matrix(0,nrow=n,ncol=m)
 
    count<-0
 
    while (delta > conv) {  
      
      ## iteration count
      count <- count + 1
      
      ## perform svd in col space
      s<-svd(X)

      ## truncate to known dimension p
      #u<-s[["u"]][,1:p]
      #d<-s[["d"]][1:p]  
      v<-s[["v"]][,1:p]
    
      ## get max likelihood estimates in row space    
     
	for (i in 1:m){
	X.maxlike.row[,i]<- v %*% solve( t(v)%*%sig.inv%*%v ) %*% t(v) %*% sig.inv %*% Xprime[,i]
      }

      ## calculate S^2
      S1<-0
      for (i in 1:m){
	S1 = S1 + t(Xprime[,i]-X.maxlike.row[,i]) %*% sig.inv %*% (Xprime[,i]-X.maxlike.row[,i])
      }

      ## perform svd using new max like estimates in row space
      s<-svd(X.maxlike.row)
    
      ## truncate
      #u<-s[["u"]][,1:p]
      #d<-s[["d"]][1:p]  
      v<-s[["v"]][,1:p]

      ## get max likelihood estimates in col space    
      for (i in 1:n){
	X.maxlike.col[,i]<- v %*% solve( t(v)%*%psi.inv%*%v ) %*% t(v) %*% psi.inv %*% X[,i]
      }
  
      ## calculate S^2
      S2<-0
      for (i in 1:n){
	S2 = S2 + t(X[,i]-X.maxlike.col[,i]) %*% psi.inv %*% (X[,i]-X.maxlike.col[,i])
      }
      
      ## set for next iteration (or return value if break)
      X<-X.maxlike.col
      Xprime<-X.maxlike.row

	## print S1,S2
	c(S1,S2)
      
      delta<-(S1-S2)/S2
  }
    
  return(list(X,Xprime,count,delta))
 }

s1<-svd(data)
s2<-svd(X)
s3<-svd(t(Xprime))

s1[[3]]
s2[[3]]
s3[[3]]

## getting same loadings up to p for orig and ML data
## fxn of variance structure? mistake in code? 

## try random data

  x1<-rnorm(10,mean=5,sd=sqrt(4))
  x2<-rnorm(10,mean=10,sd=sqrt(4))
  x3<-rnorm(10,mean=20,sd=sqrt(4))
  x4<-3*x1+2*x2
  x5<-x1-4*x2+9*x3
  x6<-x2+x3

  data.rand<-cbind(x1,x2,x3,x4,x5,x6)

  ## variance structure
  psi.rand<-diag(4,nrow=nrow(data.rand))
  sigma.rand<-diag(1,nrow=ncol(data.rand))

  mlpca.rand<-MLPCA(data.rand,psi.rand,sigma.rand,p,delta,conv)
  