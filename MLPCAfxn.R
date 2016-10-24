
#### changes on 07/31/12
## created prematrix for maximum likelihood calculation
## removed for loop used in max like calculation
## changed to one calculation each of solve(sigma) and solve(psi)
## removed for loop used in S^2 calculation

#### changes on 03/19/13
## simplified row,col notations to m,n

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
  conv<-10^-5
   
##### ALTERNATING LEAST SQUARES FOR MLPCA #####  
  
  MLPCA<-function(X,psi,sigma,p,lambda,conv){
  
    ## "column space" := original dimension of X  
    m = nrow(X)
    n = ncol(X)    

    ## transpose data
    ## call this the row space (use columns of tranpose)
    Xprime<-t(X)
   
    ## create matrices for max likelihood estimates in col and row space
    X.maxlike.col<-matrix(0,nrow = m,ncol = n)
    X.maxlike.row<-matrix(0,nrow = n,ncol = m)
 
    ## initialize iteration count 
    count<-0
    
    ## create objects for inverses of cov matrices
    ## so that solution is available immediately
    solve.sigma<-solve(sigma)
    solve.psi<-solve(psi)
   
    ## MAIN
    while (abs(lambda) > conv) {  
      
      ## iteration count
      count <- count + 1
      
      ## perform svd in col space with truncation
      s<-svd(X,nv=p)
      v<-s[["v"]]
   
      ## get max likelihood estimates in row space    
      prematrix.row<- v %*% solve( t(v)%*%solve.sigma%*%v ) %*% t(v) %*% solve.sigma
      X.maxlike.row<- as.matrix(prematrix.row) %*% as.matrix(Xprime) 

      ## calculate S^2
      diff.matrix.row<-as.matrix(Xprime-X.maxlike.row)
      S1<-sum(diag(t(diff.matrix.row)%*%solve.sigma%*%diff.matrix.row))

      ## perform svd using new max like estimates in row space with truncation
      s<-svd(X.maxlike.row,nv=p)
      v<-s[["v"]]

      ## get max likelihood estimates in col space 
      prematrix.col<- v %*% solve( t(v)%*%solve.psi%*%v ) %*% t(v) %*% solve.psi
      X.maxlike.col<- as.matrix(prematrix.col) %*% as.matrix(X)
  
      ## calculate S^2
      diff.matrix.col<-as.matrix(X-X.maxlike.col)
      S2<-sum(diag(t(diff.matrix.col)%*%solve.psi%*%diff.matrix.col))

      
      ## set for next iteration (or return value if break)
      X<-X.maxlike.col
      Xprime<-X.maxlike.row
      
      lambda<-(S2-S1)/S2
  }
    
  return(list(X,Xprime,count,lambda))
 }

system.time(data.mlpca<-MLPCA(data,psi,sigma,p,lambda,conv))




