library(rpart)		
spam.data <- read.table("http://rohan.sdsu.edu/~jjfan/sta702/spamdata.txt", 
  sep=",", na.strings="NA")
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
  labels=c("No", "Spam"))

## split data into a training sample and a test sample, using
## stratified sampling
set1<-spam.data[spam.data$spam=="Spam",]
set0<-spam.data[spam.data$spam=="Non-spam",]
dim(set1)  # 1813   58
dim(set0)  # 2788   58
> 1813*2/3
[1] 1208.667
> 2788*2/3
[1] 1858.667

# set.seed(857)  
training1<-sample(1:1813,1208)
test1<-(1:1813)[-training1]
sum((1:1813)==sort(c(training1,test1))) # 1813

training0<-sample(1:2788,1858)
test0<-(1:2788)[-training0]
sum((1:2788==sort(c(training0,test0)))) #2788

train<-rbind(set1[training1,],set0[training0,])
test<-rbind(set1[test1,], set0[test0,])
dim(train)  # 3066   58
dim(test)   # 1535   58
> 3066+1535
[1] 4601
> 1208+1858
[1] 3066

## tree growing and pruning with training data
my.control <- rpart.control(cp=0, xval=0)
fit1<- rpart(spam ~ ., data=train, method="class",
        control=my.control)

printcp(fit1)  
           CP nsplit rel error
1  0.47682119      0   1.00000
2  0.07119205      1   0.52318
3  0.06788079      3   0.38079
4  0.02235099      4   0.31291
5  0.01986755      5   0.29056
6  0.01324503      6   0.27070
7  0.00827815      8   0.24421
8  0.00745033      9   0.23593
9  0.00662252     10   0.22848
10 0.00496689     11   0.22185
11 0.00331126     14   0.20695
12 0.00248344     18   0.19371
13 0.00193157     29   0.16474
14 0.00137969     32   0.15894
15 0.00124172     35   0.15480
16 0.00082781     37   0.15232
17 0.00000000     40   0.14983

tree8<-prune(fit1,cp=(fit1$cptable[7,1]+fit1$cptable[8,1])/2)
plot(tree8,uniform=T, margin=0.2)
text(tree8,use.n=T)
pred8<-predict(tree8,newdata=test,type="class")
error8<-table(test$spam,pred8)[1,2]+table(test$spam,pred8)[2,1]

error<-rep(0,nrow(fit1$cptable))   # in the above run, nrow=17
> table(test$spam)
Non-spam     Spam 
     930      605 
error[1]<-605
for (i in 2:nrow(fit1$cptable)){
pred<-predict(prune(fit1,cp=(fit1$cptable[i-1,1]+fit1$cptable[i,1])/2),
  newdata=test,type="class")
error[i]<-table(test$spam,pred)[1,2]+table(test$spam,pred)[2,1]
}
> error
 [1] 605 328 266 218 214 207 175 167 162 156 149 138 141 138 132 134 137
> error8
[1] 167
error<-error/length(test$spam)
SE<-sqrt( error*(1-error)/length(test$spam) )
tree.n<-min( (1:nrow(fit1$cptable))[error==min(error)] ) #15

plot(1:nrow(fit1$cptable),error,type='l')
points(1:nrow(fit1$cptable),error)
abline(error[tree.n]+SE[tree.n],0)
> cbind(error,error[tree.n]+SE[tree.n])
           error          
 [1,] 0.39413681 0.0931492
 [2,] 0.21368078 0.0931492
 [3,] 0.17328990 0.0931492
 [4,] 0.14201954 0.0931492
 [5,] 0.13941368 0.0931492
 [6,] 0.13485342 0.0931492
 [7,] 0.11400651 0.0931492
 [8,] 0.10879479 0.0931492
 [9,] 0.10553746 0.0931492
[10,] 0.10162866 0.0931492
[11,] 0.09706840 0.0931492
[12,] 0.08990228 0.0931492
[13,] 0.09185668 0.0931492
[14,] 0.08990228 0.0931492
[15,] 0.08599349 0.0931492
[16,] 0.08729642 0.0931492
[17,] 0.08925081 0.0931492
## choose tree 12 by the 1 SE rule, but based on the plot 
## we may choose tree 7 if a simpler tree is desired. 
