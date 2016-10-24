train<-read.table("train.txt",header=T)
#trainNAremove<-train[-12,]
val<-read.table("val.txt",header=T)
val<-val[!is.na(val$C8),]
valpred<-val[,1:11]

require(MASS)
require(leaps)

#make response integers
A1int<-round(A1)
A2int<-round(A2)
A3int<-round(A3)
A4int<-round(A4)
A5int<-round(A5)
A6int<-round(A6)
A7int<-round(A7)

examplefit<-glm.nb(A1int~1)
y<-dnbinom(0:100,mu=exp(examplefit$coefficients[1]),size=examplefit$theta)
hist(A1int,breaks=50,freq=F)
lines(0:100,y,type="l")

#areas of non-overlap are expected to be "fixed" by addition of predictors

nbfitresults<-function(resp,newdata,observed){
  #fit glm model
  fit<-glm.nb(resp~SEASON+SIZE+VEL+C1+C2+C3+C4+C5+C6+C7+C8)
  #perform AIC selection
  select<-step(fit,direction="both",trace=F)
  #create data frame with predicted values and observed values
  validate<-as.data.frame(cbind(predict(select,newdata,type="response"),observed))
  colnames(validate)<-c("Predicted","Observed")
  #calculate sum square prediction error
  prederror<-sqrt(sum((validate$Predicted-validate$Observed)^2)/length(validate$Predicted))
  #calculate base error (estimate all zeroes)
  baseerror<-sqrt(sum((validate$Observed)^2)/length(validate$Observed))
  #plot histograms of predicted and observed
  par(mfrow=c(2,1))
  hist(validate$Predicted,breaks=50,main="Histogram of Predicted Values",xlab="Predicted Value")
  hist(validate$Observed,breaks=50,main="Histogram of Observed Values",xlab="Observed Value")
  par(mfrow=c(2,2))
  #print regression diagnostics
  plot(select)
  #fit null model (intercept only)
  nullmod<-glm.nb(resp~1)
  #get p-value for selected model
  cat("p-value from ANOVA for null model vs selected model is\n")
  print(anova(nullmod,select))
  #output
  out<-list(prederror,baseerror,anova(nullmod,select)[8],summary(select),validate)
  names(out)<-c("prederror","baseerror","p","summarySelected","predandobserved")
  return(out)
}

A1nb<-nbfitresults(A1int,valpred,val$A1)
A2nb<-nbfitresults(A2int,valpred,val$A2)
A3nb<-nbfitresults(A3int,valpred,val$A3)
A4nb<-nbfitresults(A4int,valpred,val$A4)
A5nb<-nbfitresults(A5int,valpred,val$A5)
A6nb<-nbfitresults(A6int,valpred,val$A6)
A7nb<-nbfitresults(A7int,valpred,val$A7)

A2nb$predandobserved<-A2nb$predandobserved[-46,]
A2prederror<-sqrt(sum((A2nb$predandobserved[,1]-A2nb$predandobserved[,2])^2)/nrow(A2nb$predandobserved))

A6nb$predandobserved<-A6nb$predandobserved[-c(37,46),]
A6prederror<-sqrt(sum((A6nb$predandobserved[,1]-A6nb$predandobserved[,2])^2)/nrow(A6nb$predandobserved))

RMSD<-rbind(A1nb$prederror[1],A2prederror,A3nb$prederror[1],A4nb$prederror[1],
      A5nb$prederror[1],A6prederror,A7nb$prederror[1])
rownames(RMSD)<-c("A1","A2","A3","A4","A5","A6","A7")
colnames(RMSD)<-"RMSD"

par(mfrow=c(2,1))
hist(A1nb$predandobserved[,1],breaks=50,main="",xlab="Predicted Algae Count, Species A1",xlim=c(0,80))
hist(A1nb$predandobserved[,2],breaks=50,main="",xlab="Observed Algae Count, Species A1",xlim=c(0,80))


