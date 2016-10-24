#######################
# BAGGING
#######################

library(ipred)
help(package="ipred")
?bagging    # or help(bagging)

# ==================================
# Classification: Ionosphere data
# ==================================
data(Ionosphere)
> dim(Ionosphere)
[1] 351  34
> names(Ionosphere)
 [1] "V1"    "V3"    "V4"    "V5"    "V6"    "V7"    "V8"    "V9"    "V10"  
[10] "V11"   "V12"   "V13"   "V14"   "V15"   "V16"   "V17"   "V18"   "V19"  
[19] "V20"   "V21"   "V22"   "V23"   "V24"   "V25"   "V26"   "V27"   "V28"  
[28] "V29"   "V30"   "V31"   "V32"   "V33"   "V34"   "Class"
> table(Ionosphere$V2)
  0 
351     #constant
Ionosphere$V2 <- NULL

set.seed(135)
mod <- bagging(Class ~ ., data=Ionosphere, nbagg=50, coob=T)
mod.rpart <- rpart(Class~.,data=Ionosphere)


print(mod)  #Out-of-bag estimate of misclassification error:  0.0912 

predc<-predict(mod, newdata=Ionosphere)    ##based on average of bootstrapped trees
table(predc, Ionosphere$Class)
predc  bad good
  bad  126   0 
  good   0 225 
## Note: out-of-bag estimate of error is more reliable (i.e., true error rate)
## since model is based on Ionosphere

summary(mod) # show all the details

#===============================================
# REGRESSION -- 1987 BASEBALL HITTER SALARY DATA
#===============================================

baseball<-read.table("http://www-rohan.sdsu.edu/~jjfan/sta702/bb87dat.txt", 
     header = F, col.names=c("id", "name", "bat86", "hit86", "hr86", "run86", 
  "rb86", "wlk86", "yrs", "batcr", "hitcr", "hrcr", "runcr", "rbcr", "wlkcr", 
 "leag86", "div86", "team86", "pos86", "puto86", "asst86", "err86", "salary", 
        "leag87", "team87", "logsalary"))

dim(baseball) # 263 26
set.seed(271)
fit <- bagging(logsalary~.-id - name - salary, 
   data=baseball, nbagg=50, coob=TRUE); fit
pred <- predict(fit, newdata=baseball) 


# COMPUTE R-SQUARED
n <- nrow(baseball)
sst <- (n-1)*var(baseball$logsalary)
sse <- sum((baseball$logsalary-pred)^2);
R2 <- (sst-sse)/sst
cbind(sst, sse, ssr=sst-sse, R.squared=R2)

# PLOT OF PREDICTED VS. OBSERVED
plot(baseball$logsalary, pred, ylab="predicted", xlab="observed", 
        main="Bagging Prediction for Baseball 1987")
abline(a=0, b=1, col="red")

## to assess the true prediction error, one needs to "save" some data
## from the training set that is to be used for prediction only