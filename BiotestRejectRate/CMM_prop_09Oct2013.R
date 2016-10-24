library(tolerance)

setwd("C:/users/emurphy/Dropbox/ConsultBusiness/BiotestRejectRate")

xdata = read.csv("CMM_09Oct2013.csv",header=T,sep=",")
data = xdata[!(xdata[,3] %in% c(130136,130036,130040)),]

n <- data[,6]
date <- as.Date(data[,4],"%d-%B-%Y")
crita <- apply(data[,7:15],1,sum)
critb <- apply(data[,16:17],1,sum)
maj <- apply(data[,18:21],1,sum)
min <- apply(data[,22:25],1,sum)

p.crita = crita/n
p.critb = critb/n
p.maj = maj/n
p.min = min/n

crita.utl.norm = normtol.int(p.crita,P=c(0.90,0.95,0.99),method="HE")
critb.utl.norm = normtol.int(p.critb,P=c(0.90,0.95,0.99),method="HE")
maj.utl.norm = normtol.int(p.maj,P=c(0.90,0.95,0.99),method="HE")
min.utl.norm = normtol.int(p.min,P=c(0.90,0.95,0.99),method="HE")

crita.tol = round(unlist(crita.utl.norm[5],use.names=F)*100,2)
critb.tol = round(unlist(critb.utl.norm[5],use.names=F)*100,2)
maj.tol = round(unlist(maj.utl.norm[5],use.names=F)*100,2)
min.tol = round(unlist(min.utl.norm[5],use.names=F)*100,2)


cmm.tol = rbind(crita.tol,critb.tol,maj.tol,min.tol)
rownames(cmm.tol) = c("Critical A","Critical B","Major","Minor")
colnames(cmm.tol) = c("90% Coverage","95% Coverage","99% Coverage")

write.csv(cmm.tol,"cmmtol_09Oct2013.csv")



