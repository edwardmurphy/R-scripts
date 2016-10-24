library(tolerance)

setwd("C:/users/emurphy/Dropbox/ConsultBusiness/BiotestRejectRate")

data <- read.csv("CMM_12Sep2013.csv",header=T,sep=",")

n <- data[,6]
date <- as.Date(data[,4],"%d-%B-%Y")
crit <- apply(data[,7:15],1,sum)
maj <- apply(data[,16:23],1,sum)
min <- apply(data[,24:27],1,sum)

# average values 
crit_avg = mean(crit/n)
maj_avg = mean(maj/n)
min_avg = mean(min/n)

p = c(0.90,0.95,0.99)

tolint = function(c,met){

	crit.tol = NULL
	for(i in 1:length(crit)){
		crit.tol[i] = bintol.int(round(mean(crit,0)),n[i],P=c,method=met)[[5]]
	}

	crit.true = 1-sum(crit>crit.tol)/length(crit)

	maj.tol = NULL
	for(i in 1:length(maj)){
		maj.tol[i] = bintol.int(round(mean(maj,0)),n[i],P=c,method=met)[[5]]
	}

	maj.true = 1-sum(maj>maj.tol)/length(maj)

	min.tol = NULL
	for(i in 1:length(min)){
		min.tol[i] = bintol.int(round(mean(min,0)),n[i],P=c,method=met)[[5]]
	}

	min.true = 1-sum(min>min.tol)/length(min)

	return(list(crit.tol,crit.true,maj.tol,maj.true,min.tol,min.true))
}

