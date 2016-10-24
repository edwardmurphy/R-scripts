## NOTE: Assume X: (n x r), but x_i: (r x 1)

## libraries
  library("Matrix")

## need mlpca object called
summary(data.mlpca)

## assumptions
 # observational noise matrix
  Q_eps <- sigma #same assumption as MLPCA 

## get projections from max likelihood data
  s <- svd(data.mlpca[[1]],nv=p)
  v <- s[["v"]]  #MLPCA loadings (orthonormal) 

## get emprical bayes hyperparameters
  mu_a <- as.vector(v)  # projection mean, (rp x 1)
  mu_x <- apply(data.mlpca[[1]],2,mean) # data mean, (r x 1)
  Q_x <- cov(data.mlpca[[1]]) # data covariance  (r x r)
  
  # get Q_a

    E_XTX <- Q_x + mu_x %*% t(mu_x) 
    E<-eigen(E_XTX)
    
    # create dummy list with # elements equal to number of covariance matrices
    # # elements = # of projections = # of eigenvalues of E 
    # loop through each element of list
    # within each element, loop through all projection directions
    
    n <- nrow(data.mlpca[[1]])
    cov_alpha_list <- list(NA)
    
    for (j in 1:length(E$values)){
      lam_j <- E$values[j]
      lam_k <- E$values[-j]
      alpha_j <- E$vectors[,j]
      alpha_k <- E$vectors[,-j]
      cov_alpha <- matrix(0,nrow=length(alpha_j),ncol=ncol(t(alpha_k)))
      for (i in 1:length(lam_k)){
	cov_alpha <- cov_alpha + lam_j*lam_k[i]/(lam_j-lam_k[i])^2 * alpha_j %*% t(alpha_k[,i])
      }
     cov_alpha <- 1/n * cov_alpha
     cov_alpha_list[[j]] <- cov_alpha
     }
   
   Q_a <- bdiag(cov_alpha_list[1:p])
  
    # write function in terms of unknown alpha only
    # given max like x and z from svd
    # solve for alpha
    # use new alpha to get z
    # use new alpha and z to get x
    # iterate
    
    X <- data.matrix(data)
    Z <- s[["u"]]%*%diag(s[["d"]])[,1:p]
    
    # try making into 3 functions:
    #	1. likelihood for single obs (like)
    #	2. sum of likelihoods (like_sum)
    #	3. function of 2. only to optimize ( e.g. A <- function (a) like_sum(a) )
    
    alpha_prac <- function (a) {
      sum <- 0
      alpha <- matrix(a,ncol=p,byrow=F)
      for (i in 1:nrow(X)){
	sum = sum + { t( matrix(X[i,]) - alpha%*%matrix(Z[i,]) ) %*% solve(Q_eps)%*% ( matrix(X[i,]) - alpha%*%matrix(Z[i,]) ) }
      }
    }
    
    optim(rep(0,21),alpha_prac,method="CG")
    
    
    alpha_map <- function (X, Z, Q_eps, , mu_a   ) {
      sum <- 0
      for (i in 1:nrow(X)){
	sum = sum + { t(X[i,] - matrix(a,ncol=p,byrow=F)%*%Z[i,1:p]) %*% solve(Q_eps, X[i,] - X_hat[i,]) + 
	
    
    
  
    
  








