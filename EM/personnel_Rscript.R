setwd("C:/Users/Owner/Documents/EM data")

pers<-read.table("personnel.txt",header=T)

require(MASS)

fit<-glm.nb(GownsSinceNonZero~1,data=pers)
