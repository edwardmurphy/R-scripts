###############
# RANDOM FOREST
###############

library(randomForest)
help(package="randomForest")

############
# REGRESSION
############
baseball<-read.table("http://www-rohan.sdsu.edu/~jjfan/sta702/bb87data2.txt",  
  header = F, col.names=c("bat86", "hit86", "hr86", "run86", 
  "rb86", "wlk86", "yrs", "batcr", "hitcr", "hrcr", "runcr", "rbcr", "wlkcr", 
 "leag86", "div86", "team86", "puto86", "asst86", "err86", 
        "leag87", "team87", "logsalary"))
  # data compiled by Prof. Mitch Watnik, CSU at East Bay
dim(baseball) # 263 22

# SEARCH THE BEST mtry PARAMETER, NUMBER OF VARIABLES RANDOMLY 
# SAMPLED AT EACH SPLIT 
bb.rf <- tuneRF(baseball[ , 1:21], 
   baseball[,22], ntreeTry=500, stepFactor=2, 
    improve=0.05, trace=TRUE, plot=TRUE, dobest=FALSE)
bb.rf
# the best mtry is found to be 4

# RANDOM FOREST
bb.rf <-randomForest(logsalary ~ .,
        data=baseball,  mtry=4, ntree=50, keep.forest=TRUE, importance=TRUE)
print(bb.rf) 

# LOOK AT TREE SIZES AT ALL RUNS
treesize(bb.rf, terminal=T)

# GET THE TREE AT PARTICULAR RUN
getTree(bb.rf, k=4, labelVar=TRUE) 

# LOOK AT VARIABLES ACTUALLY USED IN THE FORESTS 
varUsed(bb.rf, by.tree=FALSE)   ## this can be used as variable importance
varUsed(bb.rf, by.tree=TRUE)

# Plot the error rates or MSE of a randomForest object 
# INSPECT WHEN THE PERFORMANCE GETS STABLE
plot(bb.rf, main="MSE vs. # of boostrap Samples")

# PLOT OF VARIABLE IMPORTANCE
importance(bb.rf)
varImpPlot(bb.rf, main="Variable Importance for Baseball 1987")

################
# Classification      
################
data(iris) # data and doc available at UCI site
set.seed(71)
iris.rf <- randomForest(Species ~ ., data=iris, importance=TRUE)
print(iris.rf)

# TWO FUNCTIONS: grow AND combine
iris.rf1<-randomForest(Species~.,data=iris, ntree=50, norm.votes=F)
iris.rf2<-grow(iris.rf1,50)
print(iris.rf2)

rf1 <- randomForest(Species ~ ., iris, ntree=50, norm.votes=FALSE)
rf2 <- randomForest(Species ~ ., iris, ntree=50, norm.votes=FALSE)
rf3 <- randomForest(Species ~ ., iris, ntree=50, norm.votes=FALSE)
rf.all <- combine(rf1, rf2, rf3)
print(rf.all)

## Look at variable importance:
round(importance(iris.rf), 2)
varImpPlot(iris.rf)

## prediction
set.seed(231)
ind<-sample(2,nrow(iris),replace=T,prob=c(.8,.2)) #1/2: training/testing
iris.rf<-randomForest(Species~., data=iris[ind==1,])
iris.pred<-predict(iris.rf,iris[ind==2,])
table(observed=iris[ind==2,"Species"],predicted=iris.pred)
            predicted
observed     setosa versicolor virginica
  setosa      6      0          0       
  versicolor  0      7          0       
  virginica   0      0         13       

## compare results with a single tree
library(rpart)
my.control <- rpart.control(cp=0, xval=10)
fit1<- rpart(Species ~ ., data=iris[ind==1,], 
   method="class", control=my.control)
printcp(fit1)  
      CP nsplit rel error xerror     xstd
1 0.5375      0    1.0000 1.0625 0.064631
2 0.3875      1    0.4625 0.4625 0.063688
3 0.0000      2    0.0750 0.1125 0.036113 #best tree

besttree<-prune(fit1,cp=0.3)
plot(besttree,margin=.1)
text(besttree,use.n=T)  #tree with 3 terminal nodes
pred.tree<- predict(besttree,newdata=iris[ind==2,], type='class')
table(observed=iris[ind==2,"Species"],predicted=pred.tree)
            predicted
observed     setosa versicolor virginica
  setosa      6      0          0       
  versicolor  0      7          0       
  virginica   0      0         13         
#perfect prediction by single tree using rpart
#random forest can not do better than perfect

