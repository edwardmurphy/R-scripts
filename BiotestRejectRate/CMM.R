library(tolerance)

setwd("C:/users/emurphy/Dropbox/ConsultBusiness/BiotestRejectRate")

data <- read.csv("CMM.csv",header=T,sep=",")

n <- data[,6]
date <- as.Date(data[,4],"%d-%B-%Y")
crit <- apply(data[,7:15],1,sum)
maj <- apply(data[,16:23],1,sum)
min <- apply(data[,24:27],1,sum)

# use average values
bintol.int(round(mean(crit),0),round(mean(n),0),P=0.95)
bintol.int(round(mean(maj),0),round(mean(n),0),P=0.95)
bintol.int(round(mean(min),0),round(mean(n),0),P=0.95)

crit.tol = NULL
for(i in 1:length(crit)){
	crit.tol[i] = bintol.int(crit[i],n[i],P=0.95)[[5]]
}

