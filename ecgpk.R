library(nlme)
setwd("C:\\Users\\emurphy\\Dropbox")

data = read.csv("ecgpk.csv",header=T)
data[1:10,]
data$TAD = as.factor(data$ATPTN)
data$DAY = as.factor(data$AVISITN)

fit = lme(QTCF ~ 1 + TAD + DAY + TOTALOCA,
random = ~1 | USUBJID, data=data,na.action="na.omit")

fit1 = lme(QTCF ~ 1 + TAD + DAY + TOTALOCA, 
random = list(USUBJID = pdDiag(~TOTALOCA)),data=data,na.action="na.omit")

summary(fit1)

anova(fit1,fit)

random.effects(fit1)
intervals(fit1)
VarCorr(fit1)

predict(fit1,level=0)

fit1$apVar
fit1$modelStruct
getVarCov(fit1)
getVarCov(fit1,type="marginal")