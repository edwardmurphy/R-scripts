library(rpart)
library(randomForest) #not used
library(ipred)
library(ada)  ##response must be binary

data<-read.delim("C:\\Users\\Owner\\Documents\\Favrille\\Fav.txt",header=TRUE)
data$IgClassHCGene<-as.factor(data$IgClassHCGene)

data.lot<-data[,-c(2,8,11,24,25,26)]	#remove date index,Amp1,C2,D2,D3,fill vol,
data<-data.lot[,-1]	#remove lot

############################################################################
> colnames(data.lot)
 [1] "LotNumber"          "IgClassHCFam"       "IgClassHCGene"     
 [4] "IgClassLCFam"       "IgClassLCGene"      "NumHavestBags"     
 [7] "AvgProteinHarvBags" "C1A280"             "C4SECPurity"       
[10] "RT"                 "C6SDSRed"           "C7SDSNonRed"       
[13] "C10pH"              "HC1PropDetected"    "HC2PropDetected"   
[16] "HC3PropDetected"    "LC1PropDetected"    "LC2PropDetected"   
[19] "LC3PropDetected"    "D1BCA"              "D5ParticleSize"    
[22] "D6pH"               "D7Osmo"             "D9Endo"            
[25] "D10ResGlut"         "pI"                 "MolWtHC"           
[28] "MolWtLC"            "GlycoSitesHC"       "GlycoSitesLC"      
[31] "LysineRes"          "ExtCoeff"     
############################################################################

## tables show that nominal classes are likely insufficient to use as "responses"
## in CART, due to high # of nlevels and low sample size in some levels
table(data$IgClassHCFam)  
table(data$IgClassHCGene)
table(data$IgClassLCFam)
table(data$IgClassLCGene)

## so we need some other classification
## use MolWt as a class? create classes
quantile(2*data$MolWtHC+2*data$MolWtLC,na.rm=T)

####### kmeans cluster analysis

## use continuous responses to perform cluster analysis
## use this as classification
X<-data.lot[,6:32]
rownames(X)<-data.lot$LotNumber

#kmeans needs complete observations
X.complete<-na.omit(X)
nrow(X.complete)/nrow(X) #92% complete observations
k2<-kmeans(X.complete,2,iter.max=20,nstart=5)
k3<-kmeans(X.complete,3,iter.max=20,nstart=5)
k4<-kmeans(X.complete,4,iter.max=20,nstart=5)
k5<-kmeans(X.complete,5,iter.max=20,nstart=5)
k6<-kmeans(X.complete,6,iter.max=20,nstart=5)
k7<-kmeans(X.complete,7,iter.max=20,nstart=5)

wc2 <- sum(k2$withinss)
wc3 <- sum(k3$withinss)
wc4 <- sum(k4$withinss)
wc5 <- sum(k5$withinss)
wc6 <- sum(k6$withinss)
wc7 <- sum(k7$withinss)

plot(c(2:7),c(wc2,wc3,wc4,wc5,wc6,wc7),xlab="Number of clusters",
   ylab="Sum of squares",type="o")
lines(c(2:6),c(wc2,wc3,wc4,wc5,wc6,wc7))  ## 5 clusters

data.lot[k5$cl==5]

#hclust

hc <- hclust(dist(X.complete)^2, "cen")
plot(hc,labels=FALSE)

############ put data into process sequential groupings ###########
## this actually doesn't work well because you lose the factor
## settings which convert to dummy codes

##known from biopsy and bioinformatics
preds1<-data[,c(1:4,26:27,31)]
colnames(preds1)<-c("IgClassHCFam","IgClassHCGene","IgClassLCFam",
				"IgClassLCGene","MolWtHC","MolWtLC",
				"ExtCoeff")

##fermentation results
preds2<-data[,5:6]
colnames(preds2)<-c("NumHavestBags","AvgProteinHarvBags")

##QC purify results
preds3<-data[,7:18]
colnames(preds3)<-c("C1A20","C4SECPurity","RT","C6SDSRed","C7SDSNonRed",
				"C10pH","HC1PropDetected","HC2PropDetected",
				"HC3PropDetected","LC1PropDetected",
				"LC2PropDetected","LC3PropDetected")

##additional purify results
preds4<-data[,c(25,28:30)]
colnames(preds4)<-c("pI","GlycoSitesHC","GlycoSitesLC","LysineRes")


#FP results	
preds5<-data[,19:24]
colnames(preds5)<-c("D1BCA","D5ParticleSize","D6pH","D7Osmo",
				"D9Endo","D10ResGlut")
 
#####################################################################
getCP<-function(fit,plot.cp=FALSE,YLIM=c(0,1)){

#### Input: fit from rpart, YLIM = limits of cross validation error plot
		
#### Output: plotcp, row from cp table showing best tree from cross validation,
####         returns cp value to use as pruning value

## cptable row at which minimum xerror is found
row.xerror.min <- fit$cptable[which.min(fit$cptable[,4]),]

##calculate value.xerror.max1SE = min(xerror) + 1-SE
value.xerror.max <-  row.xerror.min[[4]] + row.xerror.min[[5]]

##save subset of cptable contains rows at which xerror is <= value.xerror.max
fit.cp.sub<-fit$cptable[fit$cptable[,4] <= value.xerror.max,]

##get cp value for pruning, cp is at max of xerror in subset
##need if statement in case only one "row" in cp table
if (length(fit.cp.sub) > 5 ) 
	cp.prune<- (fit.cp.sub[which.max(fit.cp.sub[,4]),])[[1]]
if (length(fit.cp.sub) == 5)	
	cp.prune<- fit.cp.sub[1]

if (plot.cp==TRUE) plotcp(fit)

numsplits <- fit$cptable[, 2]
trainerror <- fit$cptable[, 3]
xerror <- fit$cptable[, 4]
xstd <- fit$cptable[, 5]

plot(numsplits, trainerror, ylim=YLIM, type="l")
lines(numsplits, xerror, lty=2)
lines(numsplits, xerror-xstd, lty=3)
lines(numsplits, xerror+xstd, lty=3)
title("Cross-validation Error Estimates and Training Error")
legend("topright", c("trainerror", "xerror"), lty=c(1,2))


cat("Best Cross Validation Tree\n")
if (length(fit.cp.sub) > 5 )
	print(fit.cp.sub[which.max(fit.cp.sub[,4]),])
if (length(fit.cp.sub) == 5)	
	print(fit.cp.sub)
return(cp.prune)
}


rpartSumAndPrune<-function(fit,...  ) {
	plotcp(fit)
	rsq.rpart(fit)
	fit.prune<-prune(fit,cp=getCP(fit))
	return(fit.prune)
}


my.control <- rpart.control(xval=10, cp=0) 

##### AvgProteinHarvBags
#rpart
set.seed(3000)
fit1 <- rpart(AvgProteinHarvBags ~ IgClassHCFam+IgClassHCGene+IgClassLCFam+
		IgClassLCGene+MolWtHC+MolWtLC+ExtCoeff,data=data, 
		control=my.control,method="anova")

fit1.prune<-rpartSumAndPrune(fit1)
printcp(fit1.prune)

#bagging
set.seed(3000)
fit1.bag<-bagging(AvgProteinHarvBags ~ IgClassHCFam+IgClassHCGene+IgClassLCFam+
		IgClassLCGene+MolWtHC+MolWtLC+ExtCoeff,data=data,nbagg=50,coob=T)

fit1.pred<-predict(fit1.bag,newdata=data)
n<-nrow(data)
sst<-(n-1)*var(data$AvgProteinHarvBags,na.rm=T)
sse<-sum((data$AvgProteinHarvBags-fit1.pred)^2,na.rm=T)
R2<-(sst-sse)/sst
R2

#random forest
set.seed(3000)
data.cleaned<-na.omit(data)
fit1.tune <- tuneRF(data.cleaned[ ,c(1:4,27:28,32) ], 
   data.cleaned[,7], ntreeTry=500, stepFactor=2, 
    improve=0.05, trace=TRUE, plot=TRUE, dobest=FALSE)
fit1.tune

fit1.rf<-randomForest(AvgProteinHarvBags ~ IgClassHCFam+IgClassHCGene+IgClassLCFam+
		IgClassLCGene+MolWtHC+MolWtLC+ExtCoeff,data=data,ntree=50,
		keep.forest=TRUE,importance=TRUE,na.action=na.omit)

#####NumHavestBags
#rpart
set.seed(3000)
fit1 <- rpart(NumHavestBags ~ IgClassHCFam+IgClassHCGene+IgClassLCFam+
		IgClassLCGene+MolWtHC+MolWtLC+ExtCoeff,data=data, 
		control=my.control,method="poisson")

fit1.prune<-rpartSumAndPrune(fit1)
printcp(fit1.prune)

#bagging
set.seed(3000)
fit1.bag<-bagging(NumHavestBags ~ IgClassHCFam+IgClassHCGene+IgClassLCFam+
		IgClassLCGene+MolWtHC+MolWtLC+ExtCoeff,data=data,nbagg=50,coob=T)

fit1.pred<-predict(fit1.bag,newdata=data)
n<-nrow(data)
sst<-(n-1)*var(data$NumHavestBags,na.rm=T)
sse<-sum((data$NumHavestBags-fit1.pred)^2,na.rm=T)
R2<-(sst-sse)/sst
R2

#####RT
#rpart
set.seed(3000)
fit1 <- rpart(RT ~ IgClassHCFam+IgClassHCGene+IgClassLCFam+
		IgClassLCGene+MolWtHC+MolWtLC+ExtCoeff+NumHavestBags+AvgProteinHarvBags,data=data, 
		control=my.control,method="anova")

fit1.prune<-rpartSumAndPrune(fit1)
printcp(fit1.prune)

#bagging
set.seed(3000)
fit1.bag<-bagging(RT ~ IgClassHCFam+IgClassHCGene+IgClassLCFam+
		IgClassLCGene+MolWtHC+MolWtLC+ExtCoeff+NumHavestBags+AvgProteinHarvBags,
		data=data,nbagg=50,coob=T)

fit1.pred<-predict(fit1.bag,newdata=data)
n<-nrow(data)
sst<-(n-1)*var(data$RT,na.rm=T)
sse<-sum((data$RT-fit1.pred)^2,na.rm=T)
R2<-(sst-sse)/sst
R2

#####C7SDSNonRed
#rpart
set.seed(3000)
fit1 <- rpart(C7SDSNonRed ~ IgClassHCFam+IgClassHCGene+IgClassLCFam+
		IgClassLCGene+MolWtHC+MolWtLC+ExtCoeff+NumHavestBags+AvgProteinHarvBags,data=data, 
		control=my.control,method="anova")

fit1.prune<-rpartSumAndPrune(fit1)
printcp(fit1.prune)

#bagging
set.seed(3000)
fit1.bag<-bagging(C7SDSNonRed ~ IgClassHCFam+IgClassHCGene+IgClassLCFam+
		IgClassLCGene+MolWtHC+MolWtLC+ExtCoeff+NumHavestBags+AvgProteinHarvBags,
		data=data,nbagg=50,coob=T)

fit1.pred<-predict(fit1.bag,newdata=data)
n<-nrow(data)
sst<-(n-1)*var(data$C7SDSNonRed,na.rm=T)
sse<-sum((data$C7SDSNonRed-fit1.pred)^2,na.rm=T)
R2<-(sst-sse)/sst
R2


#####D5ParticleSize
#rpart
set.seed(3000)
fit1 <- rpart(D5ParticleSize ~ IgClassHCFam+IgClassHCGene+IgClassLCFam+
		IgClassLCGene+MolWtHC+MolWtLC+ExtCoeff+NumHavestBags+AvgProteinHarvBags+
		C1A280+C4SECPurity+RT+C6SDSRed+C7SDSNonRed+C10pH+HC1PropDetected+HC2PropDetected+
		HC3PropDetected+LC1PropDetected+LC2PropDetected+LC3PropDetected,
		data=data,control=my.control,method="anova")

fit1.prune<-rpartSumAndPrune(fit1)
printcp(fit1.prune)

#bagging
set.seed(3000)
fit1.bag<-bagging(D5ParticleSize ~ IgClassHCFam+IgClassHCGene+IgClassLCFam+
		IgClassLCGene+MolWtHC+MolWtLC+ExtCoeff+NumHavestBags+AvgProteinHarvBags+
		C1A280+C4SECPurity+RT+C6SDSRed+C7SDSNonRed+C10pH+HC1PropDetected+HC2PropDetected+
		HC3PropDetected+LC1PropDetected+LC2PropDetected+LC3PropDetected,
		data=data,nbagg=50,coob=T)

fit1.pred<-predict(fit1.bag,newdata=data)
n<-nrow(data)
sst<-(n-1)*var(data$D5ParticleSize,na.rm=T)
sse<-sum((data$D5ParticleSize-fit1.pred)^2,na.rm=T)
R2<-(sst-sse)/sst
R2

#####D10ResGlut
#rpart
set.seed(3000)
fit1 <- rpart(D10ResGlut ~ IgClassHCFam+IgClassHCGene+IgClassLCFam+
		IgClassLCGene+MolWtHC+MolWtLC+ExtCoeff+NumHavestBags+AvgProteinHarvBags+
		C1A280+C4SECPurity+RT+C6SDSRed+C7SDSNonRed+C10pH+HC1PropDetected+HC2PropDetected+
		HC3PropDetected+LC1PropDetected+LC2PropDetected+LC3PropDetected,
		data=data,control=my.control,method="anova")

fit1.prune<-rpartSumAndPrune(fit1)
printcp(fit1.prune)

#bagging
set.seed(3000)
fit1.bag<-bagging(D10ResGlut ~ IgClassHCFam+IgClassHCGene+IgClassLCFam+
		IgClassLCGene+MolWtHC+MolWtLC+ExtCoeff+NumHavestBags+AvgProteinHarvBags+
		C1A280+C4SECPurity+RT+C6SDSRed+C7SDSNonRed+C10pH+HC1PropDetected+HC2PropDetected+
		HC3PropDetected+LC1PropDetected+LC2PropDetected+LC3PropDetected,
		data=data,nbagg=50,coob=T)

fit1.pred<-predict(fit1.bag,newdata=data)
n<-nrow(data)
sst<-(n-1)*var(data$D10ResGlut,na.rm=T)
sse<-sum((data$D10ResGlut-fit1.pred)^2,na.rm=T)
R2<-(sst-sse)/sst
R2

#####GlycoSitesHC with purify results
#rpart
set.seed(3000)
fit1 <- rpart(GlycoSitesHC ~ IgClassHCFam+IgClassHCGene+IgClassLCFam+
		IgClassLCGene+MolWtHC+MolWtLC+ExtCoeff+NumHavestBags+AvgProteinHarvBags+
		C1A280+C4SECPurity+RT+C6SDSRed+C7SDSNonRed+C10pH+HC1PropDetected+HC2PropDetected+
		HC3PropDetected+LC1PropDetected+LC2PropDetected+LC3PropDetected,
		data=data,control=my.control,method="anova")

fit1.prune<-rpartSumAndPrune(fit1)
printcp(fit1.prune)

#bagging
set.seed(3000)
fit1.bag<-bagging(GlycoSitesHC ~ IgClassHCFam+IgClassHCGene+IgClassLCFam+
		IgClassLCGene+MolWtHC+MolWtLC+ExtCoeff+NumHavestBags+AvgProteinHarvBags+
		C1A280+C4SECPurity+RT+C6SDSRed+C7SDSNonRed+C10pH+HC1PropDetected+HC2PropDetected+
		HC3PropDetected+LC1PropDetected+LC2PropDetected+LC3PropDetected,
		data=data,nbagg=50,coob=T)

fit1.pred<-predict(fit1.bag,newdata=data)
n<-nrow(data)
sst<-(n-1)*var(data$GlycoSitesHC,na.rm=T)
sse<-sum((data$GlycoSitesHC-fit1.pred)^2,na.rm=T)
R2<-(sst-sse)/sst
R2

#####GlycoSitesLC
#rpart
set.seed(3000)
fit1 <- rpart(GlycoSitesLC ~ IgClassHCFam+IgClassHCGene+IgClassLCFam+
		IgClassLCGene+MolWtHC+MolWtLC+ExtCoeff+NumHavestBags+AvgProteinHarvBags+
		C1A280+C4SECPurity+RT+C6SDSRed+C7SDSNonRed+C10pH+HC1PropDetected+HC2PropDetected+
		HC3PropDetected+LC1PropDetected+LC2PropDetected+LC3PropDetected,
		data=data,control=my.control,method="anova")

fit1.prune<-rpartSumAndPrune(fit1)
printcp(fit1.prune)

#bagging
set.seed(3000)
fit1.bag<-bagging(GlycoSitesLC ~ IgClassHCFam+IgClassHCGene+IgClassLCFam+
		IgClassLCGene+MolWtHC+MolWtLC+ExtCoeff+NumHavestBags+AvgProteinHarvBags+
		C1A280+C4SECPurity+RT+C6SDSRed+C7SDSNonRed+C10pH+HC1PropDetected+HC2PropDetected+
		HC3PropDetected+LC1PropDetected+LC2PropDetected+LC3PropDetected,
		data=data,nbagg=50,coob=T)

fit1.pred<-predict(fit1.bag,newdata=data)
n<-nrow(data)
sst<-(n-1)*var(data$GlycoSitesLC,na.rm=T)
sse<-sum((data$GlycoSitesLC-fit1.pred)^2,na.rm=T)
R2<-(sst-sse)/sst
R2

