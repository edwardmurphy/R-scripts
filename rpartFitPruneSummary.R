rpartFitPruneSummary<-function(fit,YLIM=c(0,1)){

#### Input: fit from rpart, YLIM = limits of cross validation error plot
		
#### Output: plotcp, row from cp table showing best tree from cross validation,
####         returns cp value to use as pruning value

## cptable row at which minimum xerror is found
row.xerror.min <- fit$cptable[which.min(fit$cptable[,4]),]

##calculate value.xerror.max1SE = min(xerror) + 1-SE
value.xerror.max <-  row.xerror.min[[4]] + row.xerror.min[[5]]

##save subset of cptable contains rows at which xerror is <= value.xerror.max
fit.cp.sub<-fit$cptable[fit$cptable[,4] <= value.xerror.max,]

##get cp value for pruning, cp is at max of xerror in subset
cp.prune<- fit.cp.sub[1,1]

plotcp(fit)

numsplits <- fit$cptable[, 2]
trainerror <- fit$cptable[, 3]
xerror <- fit$cptable[, 4]
xstd <- fit$cptable[, 5]

plot(numsplits, trainerror, ylim=YLIM, type="l")
lines(numsplits, xerror, lty=2)
lines(numsplits, xerror-xstd, lty=3)
lines(numsplits, xerror+xstd, lty=3)
title("Cross-validation Error Estimates and Training Error")
legend("topright", c("trainerror", "xerror"), lty=c(1,2))


cat("Best Cross Validation Tree\n")
print(fit.cp.sub[which.max(fit.cp.sub[,1]),])
return(cp.prune)
}
