library(lattice)
setwd("E://Avedro//Validation")

data<-read.table("val.txt",header=T)
data$Replicate<-as.factor(data$Replicate)
data$Injection<-as.factor(data$Injection)

levels(data$Assay)<-c("Assay1","Assay2","Assay3")

xyplot(Result~Replicate|Assay,data=data,groups=data$Level,layout=c(3,1),
	aspect=1,type = "p",cex=1,pch=1,panel = function(x, ...) { 
	panel.xyplot(x,...)
	panel.abline(h=c(0.09604,0.1200,0.1441),lty=2)
	})
