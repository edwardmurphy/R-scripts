############################################################################
### Initialization
############################################################################
library(tolerance)

setwd("/home/ed/Dropbox/ConsultBusiness/BiotestRejectRate")
setwd("C:\\users\\emurphy\\Dropbox\\ConsultBusiness\\BiotestRejectRate")

#data as csv
data <- read.csv("RR_26AUG2013.csv",header=T,sep=",")

############################################################################
### function to return tolerance intervals for any subset
############################################################################

tol = function(data){
	fmrr = as.numeric(unlist(strsplit(as.character(data[,8]),split='%',fixed=TRUE)))/100
	orr = as.numeric(unlist(strsplit(as.character(data[,11]),split='%',fixed=TRUE)))/100
	N = data[,6]

	fmrr.utl.norm = normtol.int(fmrr,P=c(0.90,0.95,0.99),method="HE")
	orr.utl.norm = normtol.int(orr,P=c(0.90,0.95,0.99),method="HE")

	return(list(c(fmrr.utl.norm,orr.utl.norm)))
}

data2 = data[!(data[,4] %in% c(130136,130036)),]
data3 = data[!(data[,4] %in% c(130136,130036,130040)),]

a=tol(data)
b=tol(data2)
c=tol(data3)

ftol = round(rbind(unlist(a[[1]][5],use.names=F), unlist(b[[1]][5],use.names=F), 
	unlist(c[[1]][5],use.names=F))*100,2)

otol = round(rbind(unlist(a[[1]][10],use.names=F), unlist(b[[1]][10],use.names=F), 
	unlist(c[[1]][10],use.names=F))*100,2)

rownames(ftol) <- c("All Data (26AUG2013)","Remove lots 130136, 130036",
	"Remove lots 130136, 130036, 130040") -> rownames(otol)

colnames(ftol) <- c("90% Coverage","95% Coverage","99% Coverage") -> colnames(otol)

write.csv(ftol,"fmtol26AUG2013.csv")
write.csv(otol,"otol26AUG2013.csv")

