set.seed(3000)

library(rpart)		
spam.data <- read.table("spamdata.txt", sep=",", na.strings="NA")
dim(spam.data)  # 4601   58
names(spam.data) <- c("wfmake", "wfaddress", "wfall", "wf3d", "wfour", 
	"wfover", "wfremove", "wfinternet", "wforder", "wfmail", 
	"wfreceive", "wfwill", "wfpeople", "wfreport", "wfaddresses", 
	"wffree", "wfbusiness", "wfemail", "wfyou", "wfcredit", "wfyour", 
	"wffont", "wf000", "wfmoney", "wfhp", "wfhpl", "wfgeorge", "wf650", 
	"wflab", "wflabs", "wftelnet", "wf857", "wfdata", "wf415", "wf85", 
	"wftechnology", "wf1999", "wfparts", "wfpm", "wfdirect", "wfcs", 
	"wfmeeting", "wforiginal", "wfproject", "wfre", "wfedu", "wftable", 
	"wfconference", "cfsc", "cfpar", "cfbrack", "cfexc", "cfdollar", 
	"cfpound", "crlaverage", "crllongest", "crltotal", "spam")
table(spam.data$spam)  #0:2788; 1:  1813
spam.data$spam <- factor(spam.data$spam, levels=0:1, 
  labels=c("Non-spam", "Spam"))

### loss matrix penalizing false positive 10x higher than false negative 
lmat<-matrix(c(0,1,10,0),nrow=2,byrow=T)

### fit rpart models
fit.full<-rpart(spam~.,data=spam.data,control=rpart.control(cp=0,xval=10))
fit.loss<-rpart(spam~.,data=spam.data,control=rpart.control(cp=0,xval=10),
		parms=list(loss=lmat))


#########################################################################
EXAMPLE cp Table

1  0.4765581908      0 1.0000000 1.0000000 0.01828190
2  0.1489244346      1 0.5234418 0.5570877 0.01548619
3  0.0430226145      2 0.3745174 0.4787645 0.01463743
4  0.0308880309      4 0.2884721 0.3414231 0.01276651
5  0.0104798676      5 0.2575841 0.3022614 0.01211866
6  0.0082735797      6 0.2471042 0.2807501 0.01173553
7  0.0071704357      7 0.2388307 0.2653061 0.01144715
8  0.0052950910      8 0.2316602 0.2608935 0.01136257
9  0.0044125758     14 0.1958081 0.2388307 0.01092406
10 0.0035852179     15 0.1913955 0.2322118 0.01078711
11 0.0027578599     19 0.1770546 0.2261445 0.01065925
12 0.0025740026     22 0.1687810 0.2244898 0.01062399
13 0.0022062879     25 0.1610590 0.2333149 0.01081012
14 0.0021143593     27 0.1566464 0.2333149 0.01081012
15 0.0016547159     33 0.1439603 0.2272477 0.01068267
16 0.0011031440     36 0.1389961 0.2134584 0.01038434
17 0.0008273580     43 0.1312741 0.2195256 0.01051713
18 0.0005515720     47 0.1279647 0.2184225 0.01049317
19 0.0003677147     53 0.1246553 0.2195256 0.01051713
20 0.0000000000     62 0.1213458 0.2217319 0.01056482
##########################################################################

## use rpartFitPruneSummary function to find cp value at which to prune using 1SE rule

cp.val.full<-rpartFitPruneSummary(fit.full,c(0,1))  ### FUNCTION
cp.val.loss<-rpartFitPruneSummary(fit.loss)

##prune trees using cp value from rpartFitPruneSummary function
fit.full.prune<-prune(fit.full,cp=cp.val.full)
fit.loss.prune<-prune(fit.loss,cp=cp.val.loss)

## number of terminal nodes for full model
fit.full.prune.termnodes <- fit.full.prune$cptable[nrow(fit.full.prune$cptable),2] + 1	#23

## number of terminal nodes for loss model
fit.loss.prune.termnodes <- fit.loss.prune$cptable[nrow(fit.loss.prune$cptable),2] + 1	#35

############# plot pruned tree for full model
plot(fit.full.prune,uniform=T,margin=0.1) 
text(fit.full.prune,use.n=T)## too full so prune where number of terminal nodes is only 8

fit.full.cptable.sub <- fit.full$cptable[fit.full$cptable[,2] <= 7,]
fit.full.prune.sub <- prune(fit.full,cp=fit.full.cptable.sub[nrow(fit.full.cptable.sub),1])

pdf("full_prune.pdf")
plot(fit.full.prune.sub,uniform=T,margin=0.1)
text(fit.full.prune.sub,use.n=TRUE)
dev.off()

############# plot pruned tree for loss model
plot(fit.loss.prune,uniform=T,margin=0.1) 
text(fit.loss.prune,use.n=T)## too full so prune where number of terminal nodes is only 8

fit.loss.cptable.sub <- fit.loss$cptable[fit.loss$cptable[,2] <= 7,]
fit.loss.prune.sub <- prune(fit.loss,cp=fit.loss.cptable.sub[nrow(fit.loss.cptable.sub),1])

pdf("loss_prune.pdf")
plot(fit.loss.prune.sub,uniform=T,margin=0.1)
text(fit.loss.prune.sub,use.n=TRUE)
dev.off()


## get total error 
root.error <- min(table(spam.data$spam))/sum(table(spam.data$spam))

fit.full.total.error <- root.error * fit.full.prune$cptable[nrow(fit.full.prune$cptable),3]
fit.loss.total.error <- root.error * fit.loss.prune$cptable[nrow(fit.loss.prune$cptable),3]  ##higher than full, as expected

########## break down error into false positive and false negative

## predict
fit.full.predict <- predict(fit.full.prune,type="class")

## table with observed categorizations vs. predicted (2 x 2) 
table(spam.data$spam,fit.full.predict)

## get false positive (actual: non-spam, predict: spam) and false negative (actual: spam, predict: non-spam)
fit.full.false.positive <- table(spam.data$spam,fit.full.predict)[1,2]/sum(table(spam.data$spam,fit.full.predict))
fit.full.false.negative <- table(spam.data$spam,fit.full.predict)[2,1]/sum(table(spam.data$spam,fit.full.predict))

### do same for model with loss matrix
fit.loss.predict <- predict(fit.loss.prune,type="class")
table(spam.data$spam,fit.loss.predict)
fit.loss.false.positive <- table(spam.data$spam,fit.loss.predict)[1,2]/sum(table(spam.data$spam,fit.loss.predict))
fit.loss.false.negative <- table(spam.data$spam,fit.loss.predict)[2,1]/sum(table(spam.data$spam,fit.loss.predict))

table.error <- round(rbind(c(fit.full.false.positive,fit.full.false.negative),c(fit.loss.false.positive,fit.loss.false.negative)),3)
colnames(table.error) <- c("False Positive Rate","False Negative Rate")
rownames(table.error) <- c("1-1 Loss","10-1 Loss")









