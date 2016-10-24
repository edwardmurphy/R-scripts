
library(ada)
help(package="ada")

###########################
# EXAMPLE 1: Simulated data
###########################

n<-500;p<-10
f<-function(x,a,b,d) return( a*(x-b)^2+d )
x1 <- runif(n/2,0,4)
y1 <- f(x1,-1,2,1.7)+runif(n/2,-1,1)
x2 <- runif(n/2,2,6)
y2 <- f(x2,1,4,-1.7)+runif(n/2,-1,1)
y <- c(rep(1,n/2),rep(2,n/2))
dat <- data.frame(y=y,x1=c(x1,x2),x2=c(y1,y2), matrix(rnorm(n*8),ncol=8))
names(dat)<-c("y",paste("x",1:10,sep=""))
plot(dat$x1,dat$x2,pch=c(1:2)[y], col=c(1,8)[y], 
        xlab=names(dat)[2],ylab=names(dat)[3])
indtrain<-sample(1:n,100,replace=FALSE)
train<-dat[indtrain,]; dim(train) # 100 11
test<-dat[setdiff(1:n,indtrain),]; dim(test) # 400 11

##########################################
# DISCRETE ADABOOST under exponential loss   
##########################################

default=rpart.control()
fitdis<-ada(y~.,data=train,iter=50,loss="e",type="discrete", control=default)
print(fitdis)

# add testing data to have testing errors and testing kappa accuracies 
fitdis<-addtest(x=fitdis, test.x=test[,-1], test.y=test[,1])
summary(fitdis)
summary(fitdis,n.iter=36)

plot(fitdis,kappa=F, test=T)

pred<-predict(fitdis,train[,-1])
table(pred)

# Variable Importance
varplot(fitdis)

vip <- varplot(fitdis,plot.it=FALSE,type="scores")
round(vip,4)

###########################################
# Real AdaBoost ensemble with 4-split trees
###########################################

fitreal<-ada(y~.,data=train,iter=100,type="real",
    control=rpart.control(maxdepth=2,cp=-1,minsplit=0))
fitreal

###############################################
# Gentle AdaBoost ensemble with tree depth of 8
###############################################

fitgen<-ada(y~.,data=train,test.x=test[,-1],test.y=test[,1],iter=100,
    type="gentle",
    control=rpart.control(cp=-1,maxdepth=8))
summary(fitgen)

############################ 
# Example 2: Solubility data
############################

data(soldat)
n <- nrow(soldat)
set.seed(100)
ind <- sample(1:n)  # a permutation of 1:n
trainval <- ceiling(n*.5)
testval <- ceiling(n*.3)
train <- soldat[ind[1:trainval],]
test <- soldat[ind[(trainval+1):(trainval+testval)],]
valid <- soldat[ind[(trainval+testval+1):n],]

control <- rpart.control(cp=-1, maxdepth=14, maxcompete=1,xval=0)
gen1 <- ada(y~.,data=train,test.x=test[,-73],test.y=test[,73],
    type="gentle",control=control,iter=70)
summary(gen1)

gen2 <- addtest(gen1,valid[,-73],valid[,73])
summary(gen2) 
plot(gen2,TRUE,TRUE)
varplot(gen2)
