library(tolerance)

setwd("C:/users/emurphy/Dropbox/ConsultBusiness/BiotestRejectRate")

data <- read.csv("CMM_12Sep2013.csv",header=T,sep=",")

n <- data[,6]
date <- as.Date(data[,4],"%d-%B-%Y")
crit <- apply(data[,7:15],1,sum)
maj <- apply(data[,16:23],1,sum)
min <- apply(data[,24:27],1,sum)

p.crit = crit/n
p.maj = maj/n
p.min = min/n

crit.utl.norm = normtol.int(p.crit,P=c(0.90,0.95,0.99),method="HE")
maj.utl.norm = normtol.int(p.maj,P=c(0.90,0.95,0.99),method="HE")
min.utl.norm = normtol.int(p.min,P=c(0.90,0.95,0.99),method="HE")

crit.tol = round(unlist(crit.utl.norm[5],use.names=F)*100,2)
maj.tol = round(unlist(maj.utl.norm[5],use.names=F)*100,2)
min.tol = round(unlist(min.utl.norm[5],use.names=F)*100,2)

cmm.tol = rbind(crit.tol,maj.tol,min.tol)
rownames(cmm.tol) = c("Critical","Major","Minor")
colnames(cmm.tol) = c("90% Coverage","95% Coverage","99% Coverage")

write.csv(cmm.tol,"cmmtol_12Sep2013.csv")



