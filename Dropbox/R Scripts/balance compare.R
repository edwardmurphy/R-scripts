library(ggplot2)
setwd("E://")

options(digits=8)

HPLC<-data.frame(read.table("ALT1614.txt",header=T))
CHEM<-data.frame(read.table("ALT2006.txt",header=T))

ALL<-c(HPLC[1:nrow(HPLC),],CHEM[1:nrow(CHEM),])

summary(HPLC,digits=7)
summary(CHEM,digits=7)

par(mfrow=c(2,2))
plot(HPLC$ALT1614,ylim=c(99.9990,100.0010))
plot(CHEM$ALT2006,ylim=c(99.9990,100.0010))
plot(HPLC$ALT1614,ylim=c(99.9990,100.0010),type="l")
plot(CHEM$ALT2006,ylim=c(99.9990,100.0010),type="l")

########################################################################
##Function to plot 3 standard deviation limits and control chart limits
## based on moving range
##Inputs: data, d2=factor for MR; see Montgomery
##Outputs: plot of MR, control chart with 3sd limits and cc limits,
##		table with 3sd limits and cc limits
controlChartand3SDlims<-function(data,d2,cc)
{
  LCL.3sd<-mean(data)-3*sd(data)
  UCL.3sd<-mean(data)+3*sd(data)

  plot(data)
  abline(h=mean(data))
  abline(h=LCL.3sd,lty=2,col="blue")
  abline(h=UCL.3sd,lty=2,col="blue")
  if(cc)
  { 
    range<-rep(0,length(data)-1)
      for (i in 2:length(data))
      {
		    range[i]<-abs(data[i]-data[i-1])
	    }
    avgMR<-mean(range)
    LCL.cc<-mean(data)-3*avgMR/d2
    UCL.cc<-mean(data)+3*avgMR/d2
    abline(h=LCL.cc,lty=5,col="red")
    abline(h=UCL.cc,lty=5,col="red")
    legend(0,max(data),c("Mean","3 SD Limits","MR Limits"),lty=c(1,2,5),
	  col=c("black","blue","red"))
  }
  else
  { 
    legend(0,max(data),c("Mean","3 SD Limits"),lty=c(1,2),
    col=c("black","blue"))
  }
  
  lim.3sd<-cbind(LCL.3sd,UCL.3sd)
  colnames(lim.3sd)<-c("LCL","UCL")
  if(cc)
  {
    par(mfrow=c(1,1))
    plot(range)
    lim.cc<-cbind(LCL.cc,UCL.cc)
    lim<-rbind(lim.3sd,lim.cc)
    rownames(lim)<-c("3SD","MR")
    return(lim)
  }
  else
  {
    return(lim.3sd)
  }
    
}
#########################################################################
par(mfrow=c(2,2))
controlChartand3SDlims(HPLC$ALT1614,1.128,cc=TRUE) ##d2=1.128 b/c MR based on n=2 observations
controlChartand3SDlims(CHEM$ALT2006,1.128,cc=TRUE)

controlChartand3SDlims(ALL,1.128)

