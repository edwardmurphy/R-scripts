## read data 
data<-read.csv("U://data.csv",header=F,sep=",")
data<-read.csv("/home/ed/Dropbox/Bayesian PCA/RawData.csv",header=F,sep=",")
data<-read.csv("/home/un/m/murphy12/RawData.csv",header=F,sep=",")
data<-read.csv("C:\\Users\\emurphy\\Dropbox\\Bayesian PCA\\RawData.csv",
	header=F,sep=",")


## maturities
mat<-c(0.5,1,2,3,5,7,10)

#### exploratory

  plot(mat,data[300,],type="l",
  #ylim=c(0,.17)
  )

### sample covariance matrices
### not used yet since covariance structure is faked 
### for algorithmic learning reasons
### still to determine role of sample covariance matrix
  
  ## cov matrix of column space
  cov.col<-cov(data)

  ## cov matrix of row space
  cov.row<-cov(t(data))
  
  
#### hyperparameters of MLPCA function

  ## set size of truncation (number of components)
  p<-3

  ## COVARIANCE STRUCTURE
    ## same var for each month for any maturity (variance among values of column)
    ## value from sample variance of each column in dataset appear to all be ~0.001
    ## HOWEVER uniform variance leads to same MLPCA solution as PCA solution
    psi<-diag(0.001,nrow=nrow(data))
    
    ## show that magnitude does not alter equality of solutions
    # psi <- diag(10, nrow=nrow(data))
    
    ## add random perturbation to uniform
    cv<-0.1 #coefficient of variation
    psi.perturb.sd<-cv*mean(diag(psi))
    psi <- psi + diag(rnorm(nrow(data),mean=0, sd=psi.perturb.sd),nrow=nrow(data))
    
    ## same var for each maturity for any month (variance among values of row)
    sigma<-diag(0.00006, nrow=ncol(data))
    
    ## add random perturbation to uniform
    sigma.perturb.sd<-cv*mean(diag(sigma))
    sigma<-sigma + diag(rnorm(ncol(data),mean=0,sd = sigma.perturb.sd), nrow=ncol(data))
    
  ## convergence parameters
  lambda<-1 
  conv<-10^-2
   
##### ALTERNATING LEAST SQUARES FOR MLPCA #####  
  
  MLPCA<-function(X,psi,sigma,p,lambda,conv){
  
    ## transpose data
    ## call this the row space (use columns of tranpose)
    Xprime<-t(X)
   
    ## create matrices for max likelihood estimates in col and row space
    X.maxlike.col<-matrix(0,nrow=nrow(X),ncol=ncol(X))
    X.maxlike.row<-matrix(0,nrow=nrow(t(X)),ncol=ncol(t(X)))
 
    ## initialize iteration count and timers
    count<-0
    svdcol.time<- -1000
    maxlikecol.time<- -1000
    s1.time<- -1000
    svdrow.time<- -1000
    maxlikerow.time<- -1000
    s2.time<- -1000
    
    while (abs(lambda) > conv) {  
      
      ## iteration count
      count <- count + 1
      
      ## perform svd in col space with truncation
      svdcol.time <- rbind( svdcol.time, system.time(s<-svd(X,nv=p))[[2]] )
      u<-s[["u"]]
      d<-s[["d"]]
      v<-s[["v"]]
    
      ## get max likelihood estimates in row space    
      maxlikecol.time <- rbind( maxlikecol.time, system.time(
         for (i in 1:ncol(Xprime)){
	    X.maxlike.row[,i]<- v %*% solve( t(v)%*%solve(sigma)%*%v ) %*% t(v) %*% solve(sigma) %*% Xprime[,i]
	 })[[2]])

      ## calculate S^2
      s1.time <- rbind( s1.time, system.time(
	S1<-0
        for (i in 1:ncol(Xprime)){
	  S1 = S1 + t(Xprime[,i]-X.maxlike.row[,i]) %*% solve(sigma,Xprime[,i]-X.maxlike.row[,i])
	})[[2]])

      ## perform svd using new max like estimates in row space with truncation
      svdrow.time <- rbind( svdrow.time, system.time(s<-svd(X.maxlike.row,nv=p))[[2]] )
      u<-s[["u"]]
      d<-s[["d"]]
      v<-s[["v"]]

      ## get max likelihood estimates in col space    
      maxlikerow.time <- rbind( maxlikerow.time, system.time(
	for (i in 1:ncol(X)){
	  X.maxlike.col[,i]<- v %*% solve( t(v)%*%solve(psi)%*%v ) %*% t(v) %*% solve(psi) %*% X[,i]
	})[[2]])
  
      ## calculate S^2
      s2.time <- rbind( s2.time, system.time(
      S2<-0
	for (i in 1:ncol(X)){
	  S2 = S2 + t(X[,i]-X.maxlike.col[,i]) %*% solve(psi,X[,i]-X.maxlike.col[,i])
	})[[2]])
      
      ## set for next iteration (or return value if break)
      X<-X.maxlike.col
      Xprime<-X.maxlike.row
      
      lambda<-(S2-S1)/S2
  }
    
  return(list(X,Xprime,count,lambda,svdcol.time[-1],maxlikecol.time[-1],s1.time[-1],
	      svdrow.time[-1],maxlikerow.time[-1],s2.time[-1]))
 }

data.mlpca<-MLPCA(data,psi,sigma,p,lambda,conv)

  ## for conv = 10^-3, count = 2130
  ## 	 conv = 10^-4, count = 9422

s1<-svd(data)
s2<-svd(data.mlpca[[1]])
s3<-svd(t(data.mlpca[[2]]))

s1[[3]]
s2[[3]]
s3[[3]]

## getting same loadings up to p for orig and ML data
## fxn of variance structure? YES! when uniform

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

X<-matrix(rnorm(100,mean=10,sd=2),nrow=10)
Xerror<-X+rnorm(1)
  