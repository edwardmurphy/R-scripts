########################################################
##### project RData
save.image("ketoData.RData")
########################################################


########################################################
##### pre data stuff

source("/home/ed/Dropbox/MATH365/Project/projScript.R")

setwd("/home/ed/Dropbox/MATH365/ProjectData")
library(affy)
library(AnnotationDbi)
library(limma)
library(gplots)
########################################################



########################################################
##### read data, change labels, set control/exper groups
mydata<-ReadAffy()

#ma.labels = sampleNames(mydata)
#ma.lab = substr(ma.labels, 1, nchar(ma.labels)-4) #remove .CEL from names
#sampleNames(mydata) = ma.lab
#sampleNames(mydata)
sampleNames(mydata)<-c("NormDiet1","NormDiet2","NormDiet3","NormDiet4","NormDiet5","NormDiet6",
	"KetoDiet1","KetoDiet2","KetoDiet3","KetoDiet4","KetoDiet5")

normal.group<-mydata[,1:6]
kd.group<-mydata[,7:11]

setwd("/home/ed/Dropbox/MATH365/Project")
########################################################



########################################################
##### get images
png("image_plots%1d.png")
for (i in 1:11){
	image(mydata[,i])	
}
null<-dev.off()

##### histogram of data
png("histogram.png")
par(mfrow=c(2,1))
hist(normal.group,main="Normal Diet")
hist(kd.group,main="Ketogenic Diet")
null<-dev.off()

##### boxplot
my.colors<-rep(c("green","salmon"),each=6)
png("boxplot.png")
#boxplot(mydata,col=my.colors,las=2,horizontal=TRUE)
boxplot(mydata,col=my.colors,las=2)
null<-dev.off()

##### MA plots
png("MAplot%02d.png")
for (i in 1:11){
	for (j in (i+1):11){
		Arrayi = pm(mydata[,i])
		Arrayj = pm(mydata[,j])
		M = log2(Arrayi) - log2(Arrayj)
		A = 0.5*(log2(Arrayi) + log2(Arrayj))
		plot(A,M,pch=1,cex=0.2, main=paste("MA Plot for",sampleNames(mydata)[i],"vs",sampleNames(mydata)[j]))
		abline(h=0,col="blue")
		lines(lowess(A,M),col="red")
	}
}
null<-dev.off()

##### error rates
error.rates<-as.vector(rep(0,11))
for (i in 1:11){
	pms = pm(mydata[,i])
	mms = mm(mydata[,i])
	pm.mm = mean(pms>mms)
	pm.mm = 100*pm.mm
	error.rates[i]=100-pm.mm
}

png("MMerror.png")
plot(error.rates, type="h", main="Percentage of MMs > PMs", 
	ylab="%",xlab="Microarrays", ylim=c(0,50), col="red", lwd=5 )
grid(nx = NULL, ny = 6, col = "blue", lty = "dotted",lwd = par("lwd"),
equilogs = TRUE)
null<-dev.off()

# for latex
error.rates<-as.data.frame(error.rates)
colnames(error.rates)<-c("Error Rate")
mean.error.rate<-round(mean(error.rates[,1]),2)

##### hybridization controls

gn <- geneNames(mydata[,1])
indx.ctrls <- which(substr(gn,1,4)=="AFFX")
controls <- gn[indx.ctrls]	
hybC = controls[c(1:9,31:39)]
controls.all = log2(pm(mydata, controls))
hybC.spots = log2(pm(mydata, hybC))
arrays = sampleNames(mydata)

png("controlplot%02d.png")
for (i in 1:11){
	for (j in (i+1):11){
		plot(controls.all[,i], controls.all[,j],cex=0.7,col="blue",xlim=c(4,16),ylim=c(4,16), xlab="", ylab="")
		par(new=T)
		abline(0,1)
		par(new=T)
		plot(hybC.spots[,i], hybC.spots[,j],cex=0.7, col = "red", xlim=c(4,16),ylim=c(4,16), xlab=arrays[i], ylab=arrays[j])
		legend("topleft",legend=c("All Controls","Hybridization Controls"),
			pch=c(1,1),col=c("blue","red"),cex=c(0.7,0.7))
	}
}
null<-dev.off()
########################################################



########################################################
##### normalize replicates
data.norm <- mydata

##qnorm
qnorm.controlgroup <- normalize(normal.group,method="quantiles")
qnorm.kdgroup <- normalize(kd.group,method="quantiles")

pm(data.norm)[,1:6] <- pm(qnorm.controlgroup)
pm(data.norm)[,7:11] <- pm(qnorm.kdgroup)

mm(data.norm)[,1:6] <- mm(qnorm.controlgroup)
mm(data.norm)[,7:11] <- mm(qnorm.kdgroup)

png("rep_norm_plot_quantiles.png")
par(mfrow=c(2,1))
hist(mydata,main="Unnormalized")
hist(data.norm,main="Normalized,Quantiles")
null<-dev.off()
########################################################



########################################################
##### background correct/normalize all
##### get time for each also for speed comparison

## partial

# constant
r1<-system.time(eset.nobg.constant <- expresso(data.norm,bgcorrect.method="none",
	normalize.method="constant",pmcorrect.method="pmonly",
	summary.method="medianpolish"))

r2<-system.time(eset.bg.constant <- expresso(data.norm,bgcorrect.method="rma",
	normalize.method="constant",pmcorrect.method="pmonly",
	summary.method="medianpolish"))

# invariant set
r3<-system.time(eset.nobg.invariant <- expresso(data.norm,bgcorrect.method="none",
	normalize.method="invariantset",pmcorrect.method="pmonly",
	summary.method="medianpolish"))

r4<-system.time(eset.bg.invariant <- expresso(data.norm,bgcorrect.method="rma",
	normalize.method="invariantset",pmcorrect.method="pmonly",
	summary.method="medianpolish"))

# qspline
r5<-system.time(eset.nobg.qspline <- expresso(data.norm,bgcorrect.method="none",
	normalize.method="qspline",pmcorrect.method="pmonly",
	summary.method="medianpolish"))

r6<-system.time(eset.bg.qspline <- expresso(data.norm,bgcorrect.method="rma",
	normalize.method="qspline",pmcorrect.method="pmonly",
	summary.method="medianpolish"))

## complete

# quantile normalization
r7<-system.time(eset.nobg.quant <- expresso(data.norm,bgcorrect.method="none",
	normalize.method="quantiles",pmcorrect.method="pmonly",
	summary.method="medianpolish"))

r8<-system.time(eset.bg.quant <- expresso(data.norm,bgcorrect.method="rma",
	normalize.method="quantiles",pmcorrect.method="pmonly",
	summary.method="medianpolish"))

# loess
r9<-system.time(eset.nobg.loess <- expresso(data.norm,bgcorrect.method="none",
	normalize.method="loess",pmcorrect.method="pmonly",
	summary.method="medianpolish"))

r10<-system.time(eset.bg.loess <- expresso(data.norm,bgcorrect.method="rma",
	normalize.method="loess",pmcorrect.method="pmonly",
	summary.method="medianpolish"))

# contrasts
r11<-system.time(eset.nobg.contrast <- expresso(data.norm,bgcorrect.method="none",
	normalize.method="contrasts",pmcorrect.method="pmonly",
	summary.method="medianpolish"))

r12<-system.time(eset.bg.contrast <- expresso(data.norm,bgcorrect.method="rma",
	normalize.method="contrasts",pmcorrect.method="pmonly",
	summary.method="medianpolish"))

# no within group normalization

r13 <- system.time(eset.norepnormal.nobg.constant <- expresso(mydata,bgcorrect.method="none",
	normalize.method="constant",pmcorrect.method="pmonly",
	summary.method="medianpolish"))
########################################################



########################################################
##### access expression levels
exp.levels.nobg.constant <- exprs(eset.nobg.constant)
exp.levels.bg.constant <- exprs(eset.bg.constant)
exp.levels.nobg.invariant <- exprs(eset.nobg.invariant)
exp.levels.bg.invariant <- exprs(eset.bg.invariant)
exp.levels.nobg.qspline <- exprs(eset.nobg.qspline)
exp.levels.bg.qspline <- exprs(eset.bg.qspline)
exp.levels.nobg.quant <- exprs(eset.nobg.quant)
exp.levels.bg.quant <- exprs(eset.bg.quant)
exp.levels.nobg.loess <- exprs(eset.nobg.loess)
exp.levels.bg.loess <- exprs(eset.bg.loess)
exp.levels.nobg.contrast <- exprs(eset.nobg.contrast)
exp.levels.bg.contrast <- exprs(eset.bg.contrast)
# interesting one
exp.levels.norepnormal.nobg.constant <- exprs(eset.norepnormal.nobg.constant)
########################################################



########################################################
##### PLOTS
png("ConstantNorm.png")
par(mfrow=c(2,2))
plotDensity(exp.levels.nobg.constant,col=my.colors,main="No BG Correction, Constant Normalization")
plotDensity(exp.levels.bg.constant,col=my.colors,main="BG Correction, Constant Normalization")
boxplot(exp.levels.nobg.constant,col=my.colors,las=2)
boxplot(exp.levels.bg.constant,col=my.colors,las=2)
null<-dev.off()

png("InvariantNorm.png")
par(mfrow=c(2,2))
plotDensity(exp.levels.nobg.invariant,col=my.colors,main="No BG Correction, Invariant Set Normalization")
plotDensity(exp.levels.bg.invariant,col=my.colors,main="BG Correction, Invariant Set Normalization")
boxplot(exp.levels.nobg.invariant,col=my.colors,las=2)
boxplot(exp.levels.bg.invariant,col=my.colors,las=2)
null<-dev.off()

png("QsplineNorm.png")
par(mfrow=c(2,2))
plotDensity(exp.levels.nobg.qspline,col=my.colors,main="No BG Correction, Q Spline Normalization")
plotDensity(exp.levels.bg.qspline,col=my.colors,main="BG Correction, Q Spline Normalization")
boxplot(exp.levels.nobg.qspline,col=my.colors,las=2)
boxplot(exp.levels.bg.qspline,col=my.colors,las=2)
null<-dev.off()

png("QuantNorm.png")
par(mfrow=c(2,2))
plotDensity(exp.levels.nobg.quant,col=my.colors,main="No BG Correction, Quantile Normalization")
plotDensity(exp.levels.bg.quant,col=my.colors,main="BG Correction, Quantile Normalization")
boxplot(exp.levels.nobg.quant,col=my.colors,las=2)
boxplot(exp.levels.bg.quant,col=my.colors,las=2)
null<-dev.off()

png("LoessNorm.png")
par(mfrow=c(2,2))
plotDensity(exp.levels.nobg.loess,col=my.colors,main="No BG Correction, Loess Normalization")
plotDensity(exp.levels.bg.loess,col=my.colors,main="BG Correction, Loess Normalization")
boxplot(exp.levels.nobg.loess,col=my.colors,las=2)
boxplot(exp.levels.bg.loess,col=my.colors,las=2)
null<-dev.off()

png("ContrastNorm.png")
par(mfrow=c(2,2))
plotDensity(exp.levels.nobg.contrast,col=my.colors,main="No BG Correction, Contrasts Normalization")
plotDensity(exp.levels.bg.contrast,col=my.colors,main="BG Correction, Contrasts Normalization")
boxplot(exp.levels.nobg.contrast,col=my.colors,las=2)
boxplot(exp.levels.bg.contrast,col=my.colors,las=2)
null<-dev.off()

png("compare_repnormal.png")
par(mfrow=c(1,2))
plotDensity(exp.levels.nobg.constant,col=my.colors,main="Replicate Norm, Constant Norm")
plotDensity(exp.levels.norepnormal.nobg.constant,col=my.colors,main="No Replicate Norm, Constant Norm")
null<-dev.off()
########################################################


########################################################
##### fit limma model
design<-model.matrix(~0+factor(c(1,1,1,1,1,1,2,2,2,2,2)))
colnames(design) <- c("NormDiet", "KetoDiet")
contrast.matrix <- makeContrasts(KetoDiet - NormDiet,levels=design)
contrasts.list = c("KetoDiet vs NormDiet")

fit.nobg<-lmFit(eset.nobg.constant,design)
fit.nobg.contrast<-contrasts.fit(fit.nobg,contrast.matrix)
fit.nobg.contrast<-eBayes(fit.nobg.contrast)

fit.bg<-lmFit(eset.bg.constant,design)
fit.bg.contrast<-contrasts.fit(fit.bg,contrast.matrix)
fit.bg.contrast<-eBayes(fit.bg.contrast)
########################################################


########################################################
##### Volacano Plot
# parameters
y0=-8
y1=25
x0=-2
x1=4
i=1

M=0.27
B=0

# function
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

# plots
png("VPlot_nobg.png")
VPlotColor(fit.nobg.contrast)
null<-dev.off()

png("VPlot_bg.png")
VPlotColor(fit.bg.contrast)
null<-dev.off()
########################################################



########################################################
##### select genes based on Volcano Plot parameters
TT.nobg<-topTable(fit.nobg.contrast, coef=1, adjust = "fdr", n=length(fit.nobg.contrast))

TT.bg<-topTable(fit.bg.contrast, coef=1, adjust = "fdr", n=length(fit.bg.contrast))

selected.nobg = TT.nobg[(TT.nobg$B>B & abs(TT.nobg$logFC)>M),]
selected.bg = TT.bg[(TT.bg$B>B & abs(TT.bg$logFC)>M),]
########################################################



########################################################
##### annotate selected genes
# annota <- read.csv("RAE230A.csv", header=T, sep="\t")

"ID"->colnames(annota)[1]
genes.Info.nobg <- merge(annota,selected.nobg,by="ID")
genes.Info.bg <- merge(annota,selected.bg, by="ID")

### get particular gene info
# from annotation file
Becn1 <- annota[annota$Gene.Symbol=="Becn1",]

# from selected genes
Becn1.bg <- genes.Info.bg[genes.Info.bg$Gene.Symbol=="Becn1",]

# merge our genes with paper genes using both TT (all genes)
# and just those which are diff expressed (M >0.27)
# do MA plot for each

# how many genes match b/w annotation (15923) and paper genes (517)
annota.paper.genes <- merge(annota,paper.genes,by="Gene.Symbol")

> dim(annota.paper.genes)
[1] 320  42

# how many genes match b/w selected.nobg (1174) and paper genes (517)
genes.Info.paper.nobg <- merge(genes.Info.nobg,paper.genes,by="Gene.Symbol")

> dim(genes.Info.paper.nobg)
[1] 82 48

# how many genes match b/w selected.bg (1140) and paper genes (517)
genes.Info.paper.bg <- merge(genes.Info.bg,paper.genes,by="Gene.Symbol")

> dim(genes.Info.paper.bg)
[1] 81 48

### plots to compare our results with paper's results

# no background
plot(genes.Info.paper.nobg$Log2ratio,genes.Info.paper.nobg$logFC)

plot((genes.Info.paper.nobg$logFC+genes.Info.paper.nobg$Log2ratio)/2,genes.Info.paper.nobg$logFC-genes.Info.paper.nobg$Log2ratio,xlab="Mean",ylab="Difference")

# background
plot(genes.Info.paper.bg$Log2ratio,genes.Info.paper.bg$logFC)

plot((genes.Info.paper.bg$logFC+genes.Info.paper.bg$Log2ratio)/2,genes.Info.paper.bg$logFC-genes.Info.paper.bg$Log2ratio,xlab="Mean",ylab="Difference")

########################################################



########################################################
##### heat map
data.clus.nobg = exp.levels.nobg.constant[selected.nobg$ID,]
data.clus.bg = exp.levels.bg.constant[selected.bg$ID,]


par(oma = c(3,1,3,4))
heatmap.2(as.matrix(data.clus.nobg), col=greenred(75), scale="row", 
key=TRUE, symkey=FALSE, density.info="none", 
trace="none", cexRow=.9,cexCol=1.0, main = paste(contrasts.list[1],", Without Background Correction"))

x11()
par(oma = c(3,1,3,4))
heatmap.2(as.matrix(data.clus.bg), col=greenred(75), scale="row", 
key=TRUE, symkey=FALSE, density.info="none", 
trace="none", cexRow=.9,cexCol=1.0, main = paste(contrasts.list[1],",With Background Correction"))
########################################################





