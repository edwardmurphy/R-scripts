##The set of predictors does not appear consistent among the different types of algae (and perhaps
# no prediction for sum), so we seem to be left with a problem where we identify each algae
# type uniquely, which is beneficial for understanding how differing combinations of predictors
# lead to growth of different algae, but does not deal with fact that collected algal information 
# is a small subset of total algae community.  Still to answer- can we perform a type of discrimant
# analysis with quantitative responses? How do we classify a priori?? Can we cluster without algal 
# knowledge (so just on predictors), then see if clusters correlate to algae?

setwd("C://Users//Owner//Documents//dataexam")

Sweave("data_Sweave.nw")

require(car)
require(lattice)
require(cluster)
require(fpc)
require(leaps)
require(mi)
require(CCA)
require(MASS)
require(fields)



                                                
## ORIGINAL DATA IMPORT AND DATA SPLIT
data<-read.table("dataexam.txt",header=T)
data$C8[data$C8==10000]<-NA
# 
# rowId <- sample(nrow(data), size=127)
# train<-data[rowId,]
# val<-data[-rowId,]
# 
# write.table(train,"train.txt")
# write.table(val,"val.txt")

train<-read.table("train.txt",header=T)
#trainNAremove<-train[-12,]
val<-read.table("val.txt",header=T)
val<-val[!is.na(val$C8),]


##graphics

table(data[,1])
table(data[,2])
table(data[,3])
table(data[,1:3])

alldensityplots<-function(data){
for(i in 4:18){
  plot(density(data[,i],na.rm=T),xlab="Observed Value",main=c("Kernel Density Plot for",colnames(data)[i]))#xlab=colnames(data)[i]
}
}

pdf(file="densityplots.pdf")
par(mfrow=c(3,2))
alldensityplots(data)
dev.off()


allhistplots<-function(data){
  for (i in 4:19){
    hist(data[,i],xlab="Observed Value",main=c("Histogram for",colnames(data)[i]),breaks=40,freq=FALSE)
  }
}

pdf(file="histograms.pdf")
par(mfrow=c(3,2))
allhistplots(data)
dev.off()

logdensityplots<-function(data){
for(i in 4:19){
  plot(density(log(data[,i]),na.rm=T),xlab="Log Observed Value",main=c("Kernel Density Plot for Natural Log of",colnames(data)[i]))#xlab=colnames(data)[i]
}
}

pdf(file="logdensityplots.pdf")
par(mfrow=c(3,2))
logdensityplots(data)
dev.off()

loghistplots<-function(data){
    for (i in 4:19){
      hist(log(data[,i]),xlab="Log Observed Value",main=c("Histogram for Natural Log of",colnames(data)[i]),breaks=40,freq=FALSE)
    }
  }

pdf(file="loghistograms.pdf")
par(mfrow=c(3,2))
loghistplots(data)
dev.off()

allqqplots<-function(data){
  for (i in 4:19){
    qqnorm(data[,i],main=colnames(data)[i])
    qqline(data[,i])
    #abline(a=0,b=1)
  }
}

allqqplots(data)





#trainNAremove<-train[-12,]

## multiple imputation for training set. n=3. assumed all algae is predictor, not each individually. so 
# imputation is very much based on mean result across algaes. might get different results if each 
# algae type was treated as different dataset
imp<-mi(train)
trainimputed1<-mi.completed(imp)[[1]]
trainimputed2<-mi.completed(imp)[[2]]
trainimputed3<-mi.completed(imp)[[3]]

write.table(trainimputed1,"trainimputed1.txt")
write.table(trainimputed2,"trainimputed2.txt")
write.table(trainimputed3,"trainimputed3.txt")


## DATA EXPLORATION: graphical, cluster analysis

##lattice plot by categorical
# we don't want to do this for every response (algae type), predictor (chem) and factor(season,size,vel)
xyplot(A1~C2|SEASON,data=trainimputed1,panel = function(x, y, ...) {
  panel.xyplot(x, y, ...)
  panel.loess(x, y)
}
)

#scatter plot is very full, so need to view plots individually. use this to loop through
pdf(file="dataexamplots.pdf")
for(i in 12:18){
  
#   for(j in 1:3){
#     plot(trainimputed1[,j],trainimputed1[,i],main=colnames(trainimputed1)[j],ylab=colnames(trainimputed1)[i])
#   }
pdf(file="dataexamplots.pdf")
par(mfrow=c(4,2))
for (i in 12:12){
  for(j in 4:11){
    plot(trainNAremove[,j],trainNAremove[,i],xlab=colnames(trainNAremove)[j],ylab=colnames(trainNAremove)[i])
    lines(lowess(trainNAremove[,i]~trainNAremove[,j]))
    }
}
dev.off()

##Cluster Analysis
# current thoughts: this is interested to note, but given small sample size of cluster 2 (n=5), no further
# analysis (comparative) is possible. HOWEVER, this does give us some initial hints that differentiation among
# the algae types is possible, that different chemical concentrations leads to different algae growth, 
# not really, split is 122 v 5, not really two pops.

# pamk in fpc package determines number of clusters
chemcluster<-pamk(trainimputed1[,4:11])
algaecluster<-pamk(trainimputed1[,12:18])

chemclusterResult<-as.data.frame(chemcluster$pamobject[3])
algaeclusterResult<-as.data.frame(chemcluster$pamobject[3])

clusterResult<-cbind(chemclusterResult,algaeclusterResult)
clusterCompare<-clusterResult[,1]==clusterResult[,2]
table(clusterCompare)

trainimputed1<-cbind(trainimputed1,clusterResult)
colnames(trainimputed1)<-c(colnames(trainimputed1)[1:19],"chemcluster","algaecluster")

##clustering is the same across chem and algae, but only 5 samples in cluster 2 (out of 127). This cluster
# noted by lower C2, higher C4, Higher C5, Higher C6, Higher C7 for chem, and 
# Higher A1, Lower A2, A5, A6, A7 for algae

## plotting is not helpful, only 5 cases in one cluster

  # xyplot(A1~SEASON|chemcluster,data=train,panel = function(x, y, ...) {
  #    panel.xyplot(x, y, ...)
  #    #panel.loess(x, y)
  # }
  # )

chemK<-kmeans(trainimputed1[,4:11],2)
algaeK<-kmeans(trainimputed1[,12:18],2)

Kcompare<-chemK$cluster==algaeK$cluster
table(Kcompare)


trainimputed1<-cbind(trainimputed1,chemK$cluster)
colnames(trainimputed1)<-c(colnames(trainimputed1[1:21]),"chemKmeans")
trainimputed1$chemKmeans<-as.factor(trainimputed1$chemKmeans)

xyplot(A1~C1|chemKmeans,data=trainimputed1)

attach(train)
#attach(trainNAremove)
attach(trainimputed1)


##CCA

X<-data[,4:11] #chemical observations
Y<-log(data[,12:18]+1) #algae counts as log

correl<-matcor(X,Y)
img.matcor(correl,type=2)

#high corrs: C6/C7 (v high), C1/C8, -C2/C6, -C2/C7, C3/C7, C5/C6, C5/C7 

CC<-cc(X,Y)
plt.cc(CC,type="v",var.label=T)


Ximpute<-trainimputed1[,4:11]
Yimpute<-trainimputed1[,12:18]

cxy<-cancor(Ximpute,Yimpute)


####### REGRESSION ANALYSIS


C2sq<-C2^2
C2log<-log(C2)
C3log<-log(C3)
C4log<-log(C4)
C5log<-log(C5)
C6log<-log(C6)
C7log<-log(C7)
C8log<-log(C8+1) #need to add 1 due to presence of zero 

#negative binom GLM
attach(train)

modA1a<-glm.nb(A1~1)
modA1b<-glm.nb(A1~SEASON+SIZE+VEL+C1+C2+C3+C4+C5+C6+C7+C8)
stepA1<-step(modA1b,direction="both")
modA1c<-glm.nb(A1~VEL+C1+C3+C7)
testA1<-anova(modA1c,modA1a)

modA1d<-glm.nb(A1~SEASON+SIZE+VEL+C1+C2+C2sq+C3log+C4log+C5log+C6log+C7log+C8log)
stepA1log<-step(modA1d,direction="both")
modA1e<-glm.nb(A1~SIZE+VEL+C3log+C6log)

modA1f<-glm.nb(A1~SEASON+SIZE+VEL+C1+C2+C2sq+C3+C3log+C4+C4log+C5+C5log+C6+C6log+C7+C7log+C8+C8log)
stepA1full<-step(modA1f,direction="both")
modA1g<-glm.nb(A1~SIZE+VEL+C2+C3+C6log+C7)




modA1b<-glm.nb(A1~SEASON+SIZE+VEL+C1+C2sq+C3log+C4log+C5log+C6log+C7log+C8log)

A1log<-log(A1+1)


#lm models-
fullmodA1<-lm(A1log~SEASON+SIZE+VEL+C1+C2+C3+C4+C5+C6+C7+C8+C2sq+C2log+C3log+C4log+C5log+C6log+C7log+C8log)
  
A1select<-step(fullmodA1,direction="both",trace=F) 

fullmodASUM<-update(fullmodA1,ASUM~.)
ASUMselect<-step(fullmodASUM,direction="backward",trace=FALSE)

modA1<-lm(A1~SEASON+SIZE+VEL+C1+C2+C3+C4+C5+C6+C7+C8)
modA2<-update(modA1,A2~.)
modA3<-update(modA1,A3~.)
modA4<-update(modA1,A4~.)
modA5<-update(modA1,A5~.)
modA6<-update(modA1,A6~.)
modA7<-update(modA1,A7~.)

mod2<-update(mod1,.~.+C2sq+log(C3)-C3+log(C6)-C6+log(C7)-C7+log(C8+1)-C8)

for(j in 4:11){
  plot(train[,j],A2,xlab=colnames(train)[j])
  lines(lowess(A2~train[,j]))
}

plot(C1,ASUM)
lines(lowess(ASUM~C1),lty=2)

plot(C2,ASUM)
lines(lowess(ASUM~C2),lty=2)

plot(C3,ASUM)
lines(lowess(ASUM~C3),lty=2)

plot(C4,ASUM)
lines(lowess(ASUM~C4),lty=2)

plot(C5,ASUM)
lines(lowess(ASUM~C5),lty=2)

plot(C6,ASUM)
lines(lowess(ASUM~C6),lty=2)

plot(C7,ASUM)
lines(lowess(ASUM~C7),lty=2)

plot(C8,ASUM)
lines(lowess(ASUM~C8),lty=2)




#fit model with missing obs removed

mod1<-lm(ASUM~SEASON+SIZE+VEL+C1+C2+C3+C4+C5+C6+C7+C8)
select<-step(mod1,direction="both",trace=F)













