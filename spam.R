# This file contains some commands that might be helpful for 
# the spam problem (group project)
	
library(rpart)		
spam.data <- read.table("http://www-rohan.sdsu.edu/~jjfan/sta702/spamdata.txt", sep=",", na.strings="NA")
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

my.control <- rpart.control(cp=0, xval=10)
fit1<- rpart(spam ~ ., data=spam.data, method="class",
        control=my.control)

printcp(fit1)  # 0.21236 +0.010360= 0.22272, tree 15 with nsplit=33
plotcp(fit1)  # this is a typical example where xerror first drops rapidly
              # and then has a long flat valley

#tree 15
fit1pruned <- prune(fit1, cp=.0017)
plot(fit1pruned, uniform=T)
text(fit1pruned, use.n=T)  #too crowded to see

#tree 7 with 8 terminal nodes
fit1pr2 <- prune(fit1, cp=.008)
plot(fit1pr2, uniform=T, margin=0.1)
text(fit1pr2, use.n=T)  #much better

#make predictions using the optimal tree
fitted1 <- predict(fit1pruned, type="class")
> table(spam,fitted1)
    fitted1
spam Non-spam Spam
   0 2671      117
   1  144     1669

### try penalizing false positive errors more heavily
loss2 <- matrix(c(0,2,1,0), ncol=2, byrow=F)
fit2 <- rpart(spam ~ ., data=spam.data, method="class",
	control=my.control, parms=list(loss=loss2))

printcp(fit2)  # 0.27192 +0.013956= 0.285876 #tree 13 with nsplit=23

plotcp(fit2)  # again, xerror first drops rapidly
              # and then has a long flat valley

#tree 13
fit2pruned <- prune(fit2, cp=.0034)
plot(fit2pruned, uniform=T)
text(fit2pruned, use.n=T)  #too crowded to see

#tree 7 with 8 terminal nodes
fit2pr2 <- prune(fit2, cp=.012)
plot(fit2pr2, uniform=T, margin=0.1)
text(fit2pr2, use.n=T)  #much better

#make predictions using the optimal tree
fitted2 <- predict(fit2pruned, type="class")
> table(spam,fitted2)
    fitted2
spam Non-spam Spam
   0 2730       58
   1  269     1544
