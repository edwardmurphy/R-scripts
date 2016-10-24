#function
signTestUpperAlt <- function(x,theta0,p,level){
	n <- length(x)	
	sum_psi <- sum(x>theta0)
	alph <- sum(dbinom(sum_psi:n,n,p))
	if(alph>level) paste(cat("The null hypothesis is NOT rejected,alpha=",alph,"\n"))
	else paste(cat("The null hypothesis is rejected,alpha=",alph,"\n"))
	return(alph)
}

#data
x<-c(36,31,30,27,20,33,27,18,19,28)

#test
pr2.2.3 <- signTestUpperAlt(x,22,0.5,0.0107)

