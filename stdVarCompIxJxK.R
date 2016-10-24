library(lattice)

#############################################################################
##### Script to perform ANOVA and lme modeling (id.) for qual/validation data
##### from DOE using i assays/j levels/k injections
#############################################################################

#############################################################################
#### Data must be in varCompAnalysisR.txt file in working dir with 
#### the following headers (col order can be arbitrary, as cols are called 
#### by name):
#### column 1: Assay
#### column 2: Level (Low, Mid, High)
#### column 3: Replicate
#### column 4: Injection
#### column 5: Result

#### Missing data should be filled in prior to reading with 0 (zero)
#### If 0 (zero) is a legitimate value, change script to read dummy value
#### chosen 

#### Need expected results in expected.txt file in working dir with
#### expected results for each level in a column vector
#############################################################################

stdVarComp <- function ( ) {

	# Read in data, clean up, factorize levels appropriately
	data<-read.table("varCompAnalysisR.txt",header=T)
	data[data==0]<-NA   #missing data

	data$Assay<-as.factor(data$Assay)
	levels(data$Assay)<-c("Assay1","Assay2","Assay3")

	data$Level<-as.factor(data$Level)
	levels(data$Level)<-c("Low","Mid","High")

	data$Replicate<-as.factor(data$Replicate)
	data$Injection<-as.factor(data$Injection)

	# Read in expected results
	expect <- as.vector(read.table("expected.txt"))

	# plot results using lattice plot
	xyplot(Result~Replicate|Assay,data=data,groups=data$Level,
		layout=c(nlevels(data$Assay,1),aspect=1,type = "p",cex=1,#pch=c(0,1,2),
		panel = function(x, ...) { 
			panel.xyplot(x,...)
			panel.abline(h=c(LowLevel,MidLevel,HighLevel),lty=2)
			# how to generalize adding lines???   
		}
		)

	# create recovery tables for each replicate
	# view recovery to nearest 0.1%

	# create array to hold recovery results by replicate
	# first dimension (row) holds assay result
	# second dimension (col) holds replicate result
	# third dimension (holds level result)
	# so RecReplicate[,,1] is all recovery results for level 1
	# create 2 arrays- one with unrounded results, the other with rounded results (for display)

	RecReplicate <- array(dim=c(nlevels(data$Assay),nlevels(data$Replicate),nlevels(data$Level)))
	for (i in 1:nlevels(data$Level)){
		for (j in 1:nlevels(data$Assay)){

			RecReplicate[j,,i] <- with (subset(data,Level==levels(data$Level)[[i]]&Assay==levels(data$Assay)[[j]]),
					tapply(Result, Replicate, mean)/expect[i,1] * 100)

		}
	}

	RecReplicateRound <- round(RecReplicate,1)

	# create recovery table for each assay
	# row contains assay mean
	# col contains level
	# so first row contains the mean results for assay 1 at each level	

	RecAssay <- apply(RecReplicate,c(1,3),mean)
	RecAssayRound <- round(RecAssay,1)

	# create recovery table for each level

	RecLevel <- apply(RecAssay,2,mean)
	RecLevelRound <- round(RecLevel,1)

	# perform ANOVA

	for (i in 1:nlevels(data$Level)) {

		aov.mod<-aov(Result~1+Error(Assay/Replicate),data=subset(data,data$Level=="Low"))
		
		str(summary(aov.mod))
		



	



