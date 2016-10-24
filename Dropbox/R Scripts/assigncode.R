### R code from vignette source 'assign.Rnw'

###################################################
### code chunk number 1: assign.Rnw:56-71
###################################################
setwd("/home/ed/Dropbox/MATH365/Project")
load("ketoData.RData")
setwd("/home/ed/Dropbox/MATH365")
design<-model.matrix(~0+factor(c(1,1,1,1,1,1,2,2,2,2,2)))
colnames(design) <- c("NormDiet", "KetoDiet")
contrast.matrix <- makeContrasts(NormDiet-KetoDiet,levels=design)
contrasts.list = c("NormDiet vs KetoDiet")

fit.nobg<-lmFit(eset.nobg.constant,design)
fit.nobg.contrast<-contrasts.fit(fit.nobg,contrast.matrix)
fit.nobg.contrast<-eBayes(fit.nobg.contrast)

fit.bg<-lmFit(eset.bg.constant,design)
fit.bg.contrast<-contrasts.fit(fit.bg,contrast.matrix)
fit.bg.contrast<-eBayes(fit.bg.contrast)


###################################################
### code chunk number 2: assign.Rnw:74-112
###################################################
y0=-8
y1=25
x0=-2
x1=4

M=1
B=0
i=1
VPlotColor <- function(fit2){
		i=1
    volcanoplot(fit2, col="blue", ylim=c(y0,y1), xlim=c(x0,x1), cex.lab=1.3,main="")
    par(new=T)
    abline(v=-M, col="brown")
    par(new=T)
    abline(v=M, col="brown")
    par(new=T)
    abline(h=B, col="black")
    par(new=T)
    ind1 = abs(fit2$coef[,i])>M
    ind2 = fit2$lods[,i] >B
    ind3 = (fit2$coef[,i]>M & fit2$lods[,i]>B)
    ind4 = (fit2$coef[,i]< -M & fit2$lods[,i]>B)
    x = as.matrix(fit2$coef[ind1,i])
    y = as.matrix(fit2$lods[ind1,i])
    plot(x, y, col="magenta",ylim=c(y0,y1), xlim=c(x0,x1),main="",xlab="",ylab="",pch = "*",cex.lab=1.2)
    x = as.matrix(fit2$coef[ind2,i])
    y = as.matrix(fit2$lods[ind2,i])
    par(new=T)
    plot(x, y, col="orange",  ylim=c(y0,y1), xlim=c(x0,x1), main="", pch = 19,xlab="", ylab="",cex.lab=1.2)
    x = as.matrix(fit2$coef[ind3,i])
    y = as.matrix(fit2$lods[ind3,i])
    par(new=T)
    plot(x, y, col="red",  ylim=c(y0,y1), xlim=c(x0,x1),main="", pch = 19,xlab="", ylab="",cex.lab=1.2)
    x = as.matrix(fit2$coef[ind4,i])
    y = as.matrix(fit2$lods[ind4,i])
    par(new=T)
    plot(x, y, col="darkgreen",  ylim=c(y0,y1), xlim=c(x0,x1),main="", pch = 19,xlab="", ylab="",cex.lab=1.2)
 }


###################################################
### code chunk number 3: assign.Rnw:115-118
###################################################
png("VPlot_nobg.png")
VPlotColor(fit.nobg.contrast)
null<-dev.off()


###################################################
### code chunk number 4: assign.Rnw:121-124
###################################################
png("VPlot_bg.png")
VPlotColor(fit.bg.contrast)
null<-dev.off()


###################################################
### code chunk number 5: assign.Rnw:155-156
###################################################
xtable(topTable(fit.nobg.contrast, coef=1, adjust = "BH", n=40),digits=c(3,3,3,3,3,-3,-3,3), caption= "Summary statistics for linear fit of expression set with contrast normalization and no background correction: 40 most significant genes by B statistic", label="tab:nobg")


###################################################
### code chunk number 6: assign.Rnw:159-160
###################################################
xtable(topTable(fit.bg.contrast, coef=1, adjust = "BH", n=40),digits=c(3,3,3,3,3,-3,-3,3), caption= "Summary statistics for linear fit of expression set with contrast normalization and RMA background correction: 40 most significant genes by B statistic", label="tab:bg")


