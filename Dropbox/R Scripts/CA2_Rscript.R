twobinoms<-function(a1, b1, y1, n1, a2, b2, y2, n2, n.size,graph)
{
  p1<-rbeta(n.size,a1+y1, b1+n1-y1)
  p2<-rbeta(n.size,a2+y2, b2+n2-y2)
  dp1p2<-p1-p2
  Pp1bigger<-sum(dp1p2>0)/n.size
  varP<-(Pp1bigger)*(1-Pp1bigger)/n.size
  meand<-mean(dp1p2)
  medd<-median(dp1p2)
  CI<-quantile(dp1p2, c(.025, .975))
  if(graph)
  {
    plot(density(dp1p2), main="Kernel Density Estimate of p1-p2", xlab="", ylab="")
  }
  output<-c(meand, medd, Pp1bigger, CI, varP)
  names(output)<-c("mean diff", "med diff", "diff>0", ".025", ".975","Var")
  return(output)
}  

postStats<-twobinoms(2,1,8,10,1,1,13,20,500,graph=T)
priorProb<-sum((rbeta(500,2,1)-rbeta(500,1,1))>0)/500
BF<-(postStats[[3]]/(1-postStats[[3]]))/(priorProb/(1-priorProb))
BF


twobinomsrepeat<-function(iter,a1,b1,y1, n1, a2, b2, y2, n2, n.size,...)
{
  iterDiff<-rep(0,times=iter)
  iterVar<-rep(0,times=iter)
  for(i in 1:iter)
  ####figure out way to do this without for loop
    ##create a matrix where the output from twobinoms goes to and apply 
    ##this across each column? 
  {
    iterRun<-twobinoms(a1, b1, y1, n1, a2, b2, y2, n2, n.size,graph=FALSE,...)
    iterDiff[i]<-iterRun[[3]]
    iterVar[i]<-iterRun[[6]]
  }
  plot(density(iterDiff))
  return(cbind(iterDiff,iterVar))
}

diffrepeat<-twobinomsrepeat(1000,2,1,8,10,1,1,13,20,500)




