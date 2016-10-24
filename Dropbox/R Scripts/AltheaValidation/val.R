library(lattice)
setwd("E://Santen//Validation")

data<-read.table("val.txt",header=T)
data$Replicate<-as.factor(data$Replicate)
data$Injection<-as.factor(data$Injection)

levels(data$Assay)<-c("Assay1","Assay2","Assay3")

xyplot(Amount~Replicate|Assay,data=data,groups=data$Level,layout=c(3,1),
	aspect=1,type = "p",cex=1,pch=c(0,1,2),col=c("black","black","black"),
	panel = function(x, ...) { 
	panel.xyplot(x,...)
	#panel.abline(h=c(0.09604,0.1200,0.1441),lty=2)
	})

update(trellis.last.object(),
	key = list(text = list(c("80% Level", "100% Level", "120% Level")),
		space="top",
		points = list(cex = 1, pch = c(0,1,2),
            type = "p"), title="Legend"))

