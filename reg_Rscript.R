prederror<-numeric(ncol(OBS))
for(i in 1:ncol(OBS)){
  diffsq<-(PR[,i]-OBS[,i])^2
  sumsq<-sum(diffsq)
  prederror[i]<-sqrt(sumsq/length(OBS[,i]))
}



require(leaps)



attach(trainNAremove)

C1sq<-C1^2
C2sq<-C2^2
C3sq<-C3^2
C4sq<-C4^2
C5sq<-C5^2
C6sq<-C6^2
C7sq<-C7^2
C8sq<-C8^2
C2log<-log(C2)
C3log<-log(C3)
C4log<-log(C4)
C5log<-log(C5)
C6log<-log(C6)
C7log<-log(C7)
C8log<-log(C8+1)

A1logPlus1<-log(A1+1)
A2logPlus1<-log(A2+1)
A3logPlus1<-log(A3+1)
A4logPlus1<-log(A4+1)
A5logPlus1<-log(A5+1)
A6logPlus1<-log(A6+1)
A7logPlus1<-log(A7+1)

resplog<-as.matrix(cbind(A1logPlus1,A2logPlus1,A3logPlus1,A4logPlus1,A5logPlus1,
    A6logPlus1,A7logPlus1))
    
    
fitAll<-function(resp){
fitnull<-lm(resp~1)
#all predictors with transforms includes + original form
fitall<-lm(resp~SEASON+SIZE+VEL+C1+C1log+C1^2+C2+C2log+C2sq+C3+C3log+C3sq+
        C4+C4log+C4sq+C5+C5log+C5sq+C6+C6log+C6sq+
        C7+C7log+C7sq+C8+C8log+C8sq)
#model selection
selectall<-step(fitall,direction="both",trace=F)
#print(anova(fitnull,selectall))
cat("Prediction Error\n")
print(sqrt(sum((selectall$fitted.values-resp)^2)/length(resp)))
return(selectall)
}

fitSub<-function(resp){
fitnull<-lm(resp~1)
#subset of transformed predictors based on data exploration        
fitsub<-lm(resp~SEASON+SIZE+VEL+C1+C2+C2sq+C3log+C4log+C5log+C6log+C7log+C8log)
#model selection
selectsub<-step(fitsub,direction="both",trace=F)
print(anova(fitnull,selectsub))
cat("Prediction Error\n")
print(sqrt(sum((selectsub$fitted.values-resp)^2)/length(resp)))
return(selectsub)
}

# A1all<-fitAll(A1logPlus1)
# A1sub<-fitSub(A1logPlus1)
# 
# A2all<-fitAll(A2logPlus1)
# A2sub<-fitSub(A2logPlus1)
# 
# A3all<-fitAll(A3logPlus1)
# A3sub<-fitSub(A3logPlus1)
# 
# A4all<-fitAll(A4logPlus1)
# A4sub<-fitSub(A4logPlus1)
# 
# A5all<-fitAll(A5logPlus1)
# A5sub<-fitSub(A5logPlus1)
# 
# A6all<-fitAll(A6logPlus1)
# A6sub<-fitSub(A6logPlus1)
# 
# A7all<-fitAll(A7logPlus1)
# A7sub<-fitSub(A7logPlus1)

a1all<-fitAll(A1)
a1sub<-fitSub(A1)

a2all<-fitAll(A2)
a2sub<-fitSub(A2)

a3all<-fitAll(A3)
a3sub<-fitSub(A3)

a4all<-fitAll(A4)
a4sub<-fitSub(A4)

a5all<-fitAll(A5)
a5sub<-fitSub(A5)

a6all<-fitAll(A6)
a6sub<-fitSub(A6)

a7all<-fitAll(A7)
a7sub<-fitSub(A7)

# A1nb<-glm.nb(A1~SEASON+SIZE+VEL+C1+C2+C3+C4+C5+C6+C7+C8)
# A1nbstep<-step(A1nb,direction="both",trace=F)
# sqrt(sum((predict(A1nbstep,type="response")-A1)^2)/length(A1))
# 
# A2nb<-glm.nb(A2~SEASON+SIZE+VEL+C1+C2+C3+C4+C5+C6+C7+C8)
# A2nbstep<-step(A2nb,direction="both",trace=F)
# sqrt(sum((predict(A2nbstep,type="response")-A2)^2)/length(A2))

valpred<-val[,1:11]
valpred$C1sq<-(valpred$C1)^2
valpred$C2sq<-(valpred$C2)^2
valpred$C3sq<-(valpred$C3)^2
valpred$C4sq<-(valpred$C4)^2
valpred$C5sq<-(valpred$C5)^2
valpred$C6sq<-(valpred$C6)^2
valpred$C7sq<-(valpred$C7)^2
valpred$C8sq<-(valpred$C8)^2

valpred$C1log<-log(valpred$C1)
valpred$C2log<-log(valpred$C2)
valpred$C3log<-log(valpred$C3)
valpred$C4log<-log(valpred$C4)
valpred$C5log<-log(valpred$C5)
valpred$C6log<-log(valpred$C6)
valpred$C7log<-log(valpred$C7)
valpred$C8log<-log(valpred$C8+1)

a1val<-as.data.frame(cbind(predict(a1all,valpred),val[,12]))
colnames(a1val)<-c("predict","observe")
a1val$predict[a1val$predict<0]<-0
a1error<-sqrt(sum((a1val$predict-a1val$observe)^2)/length(a1val$predict))

a2val<-as.data.frame(cbind(predict(a2all,valpred),val[,13]))
colnames(a2val)<-c("predict","observe")
a2val$predict[a2val$predict<0]<-0
a2error<-sqrt(sum((a2val$predict-a2val$observe)^2)/length(a2val$predict))
