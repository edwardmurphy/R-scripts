###Calculate power for hypergeometric sampling
###Want beta ~= 0.05 at 0 observations (failures)

a<-rhyper(100000,0.1*1300,0.9*1300,30)
cumsum(table(a)/length(a))

a<-rhyper(100000,0.03*1300,0.97*1300,92)
cumsum(table(a)/length(a))

a<-rhyper(100000,0.05*1300,0.95*1300,56
cumsum(table(a)/length(a))
