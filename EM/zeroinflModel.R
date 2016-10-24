###########  understand zerioinfl output ############
r = rzinbinom(10000,mu=125,size=0.2,zprob=0.2)
sex = c(rep("MALE",5000),rep("FEMALE",5000))
d = as.data.frame(cbind(r,sex))
for(i in 1:length(r)){d$resn[i] = rzinbinom(1,mu=125,size=0.2,zprob=0.2)}
test = zeroinfl(formula = resn~sex|sex, dist="negbin",data=d,EM=T)
summary(test)
exp(coef(test)[1]) #mu
plogis(coef(test)[2]) #zprob
newdata = expand.grid(list(sex="MALE"))
phat = predprob(test,newdata)
p = as.data.frame(cumsum(phat))
male = subset(d,d$sex=="MALE")
cutoff=0.6

as.numeric(rownames(subset(p,p>=cutoff,drop=FALSE)))[1] - 1
quantile(male$resn,prob=cutoff)

#####################################################

########## libraries, working directory ###############
library(pscl)
library(lattice)
library(plyr)
library(emdbook)
setwd("C:\\Users\\emurphy\\Dropbox\\ConsultBusiness\\EM")
#######################################################

############ read data, change factors and summarize ###############
data12 = read.table("TP2012.txt",header=T)
data12$YEAR = "2012"
data13 = read.table("TP2013.txt",header=T)
data13$YEAR = "2013"

## incorrect MEAS (room) for some observations
data12$MEAS[data12$MEAS=="527DL"] = "TP"

data12$PARAM = as.factor(data12$PARAM)
data13$PARAM = as.factor(data13$PARAM)
#data12$PARAM = revalue(data12$PARAM,c("0.5"="low","5"="high"))

# one outlier result in 2013 data
outly = data13[data13$RES==max(data13$RES,na.rm=T),]
data13 = data13[data13$RES<max(data13$RES,na.rm=T),]

summary(data12)
summary(data13)
#levels(data12$ROOM)
#levels(data12$SITE)
#levels(data12$RES)

data.all = rbind(data12,data13)

class.list = c("ISO5","ISO6","ISO7","ISO8")
param.list = c("0.5","5")
######################################################################

################## proportion of zeroes in each class/param combo ###########################
prop.zero = function(data,class,param){
	sub = subset(data,CLASS==class & PARAM==param)
	sub.zero = subset(sub,RES==0)
	p = length(sub.zero$RES)/length(sub$RES)
	return(p)
}

propzero = function(data){
	props = matrix(nrow=length(class.list),ncol=length(param.list))
	for (i in 1:length(class.list)){
		for (j in 1:length(param.list)){
			props[i,j]=prop.zero(data,class.list[i],param.list[j])
		}

	}
	colnames(props) = param.list
	rownames(props) = class.list
	return(props)
}

propzero12 = propzero(data12)
propzero13 = propzero(data13)
propzeroall = propzero(data.all)

propzero12
propzero13
propzeroall
#ISO6 changed dramatically
###############################################################################################

###################### function to get quantiles #####################################
quants = function(data,j){
	q = matrix(nrow=length(class.list),ncol=length(param.list))
	for (i in 1:length(class.list)){
		sub = subset(data,CLASS==class.list[i] & PARAM==param.list[j])
		sub.zero = subset(sub,RES!="NA")
		q[i,]=quantile(sub.zero$RES,prob=c(0.95,0.99))
	}
	rownames(q)=class.list
	colnames(q)=c("95%","99%")
	return(q)
}

quants12 = list(quants(data12,1),quants(data12,2))
names(quants12) = c("0.5 micron","5 micron")
quants13 = list(quants(data13,1),quants(data13,2))
names(quants13) = c("0.5 micron","5 micron")
quantsall = list(quants(data.all,1),quants(data.all,2))

quants12
quants13
quantsall
#not dramatically different, ISO5 0.5um lower in 2013
########################################################################################


############################## PLOTS ###########################################
xyplot(RES~CLASS|PARAM,data=data12)
xyplot(RES~CLASS|PARAM,data=data13)
xyplot(RES~CLASS|YEAR+PARAM,data=data.all)
#################################################################################

########################### MODELS #############################################

mod12.zinb = zeroinfl(formula = RES~CLASS+PARAM|CLASS+PARAM, data=data12, dist="negbin",EM=T)
mod13.zinb = zeroinfl(formula = RES~CLASS+PARAM|CLASS+PARAM, data=data13, dist="negbin",EM=T)
modall.zinb = zeroinfl(formula = RES~CLASS+PARAM|CLASS+PARAM, data=data.all, dist="negbin",EM=T)

mod12.nb = glm.nb(RES ~ CLASS+PARAM, data = data12)
mod13.nb = glm.nb(RES ~ CLASS+PARAM, data = data13)
modall.nb = glm.nb(RES ~ CLASS+PARAM, data = data.all)

##### model comparisons
vuong(mod12.zinb, mod12.nb)
vuong(mod13.zinb, mod13.nb)
vuong(modall.zinb, modall.nb) 

# is year significant?
modyear.zinb = zeroinfl(formula = RES~CLASS+PARAM+YEAR|CLASS+PARAM+YEAR, data=data.all, dist="negbin",EM=T)
summary(modyear.zinb)

# add room 
modsite.zinb = zeroinfl(formula = RES~CLASS+PARAM+SITE|CLASS+PARAM+SITE, data=data.all, dist="negbin",EM=T)

#mod4 = zeroinfl(formula = RES~CLASS|1, data=data12, dist="poisson",EM=T)
#mod5 = glm(RES~CLASS, family=poisson,data=data12)
#vuong(mod4, mod5)
#################################################################################

################# PREDICTED PROBABILITIES #########################################################

##### CALCULATION FOR ISO5, PARAM = 0.5, YEAR = 2012
zprob = plogis(-1.82450) 
mu = exp(3.20388) 
size = 0.3135 
iso5=rzinbinom(10000,mu=mu,size=size,zprob=zprob)

preds = function(model,class,param,cutoff){
	predcount = matrix(nrow=length(class),ncol=length(param))
	probzero = matrix(nrow=length(class),ncol=length(param))
	for(i in 1:length(class)){
		for(j in 1:length(param)){
			newdata = expand.grid(list(CLASS=class[i],PARAM=param[j]))
			phat <- predprob(model, newdata = newdata)
			probzero[i,j] = phat[1,1]
			p = as.data.frame(cumsum(phat))
			predcount[i,j] = as.numeric(rownames(subset(p,p>=cutoff,drop=FALSE)))[1] - 1
		}
	}
	cat("Probability = ",cutoff,"\n")
	colnames(predcount) = param.list
	rownames(predcount) = class.list
	#return(list(predcount,probzero))
	return(predcount)
}

al2012=preds(mod12.zinb,class.list,param.list,0.95)
ac2012=cbind(al2012,preds(mod12.zinb,class.list,param.list,0.99))
a
al2013=preds(mod13.zinb,class.list,param.list,0.95)
ac2013=cbind(al2013,preds(mod13.zinb,class.list,param.list,0.99))

alall=preds(modall.zinb,class.list,param.list,0.95)
acall=cbind(alall,preds(modall.zinb,class.list,param.list,0.99))

write.table(ac2012,file="C:\\Users\\emurphy\\Dropbox\\ConsultBusiness\\EM\\tplimits2012.csv")
write.table(ac2013,file="C:\\Users\\emurphy\\Dropbox\\ConsultBusiness\\EM\\tplimits2013.csv")
write.table(acall,file="C:\\Users\\emurphy\\Dropbox\\ConsultBusiness\\EM\\tplimitsoverall.csv")

############ PLOTS OF ACTUAL AND PREDICTED ##########

histoverlay = function(data,class,param,model,breaksize,upplim){
	sub = subset(data,CLASS==class & PARAM==param)
	hist(sub$RES,breaks=breaksize,freq=F,xlim=c(0,upplim),xlab=paste("Class = ",class,"Size = ",param))
	newdata = expand.grid(list(CLASS=class,PARAM=param))
	phat <- predprob(model, newdata = newdata)
	m = length(phat)
	count = seq(0,m-1,by=1)
	phat2 = data.frame(cbind(as.vector(phat),count))
	lines(phat2$count,phat2$V1,type="l")
}

histoverlay(data12,"ISO5","0.5",mod12.zinb,5000,100)
histoverlay(data12,"ISO5","5",mod12.zinb,1000,20)



############################ RESAMPLING METHOD WITH MODEL PARAMS TO GET QUANTILES ################################
#### does not work!!! 

nbprobs = function(model,j){
	q.mod = matrix(nrow=length(class.list),ncol=length(param.list))
	coeff = model$coefficients
	theta = model$theta
	for (i in 1:length(class.list)){
		if (j==1) {
			if (i==1) {
				zprob = plogis(coeff$zero[1])
				mu = exp(coeff$count[1])
				r = rzinbinom(100000,mu=mu,size=theta,zprob=zprob)
				q.mod[i,] = quantile(r,prob=c(0.95,0.99))
			}
			else {
				zprob = plogis(coeff$zero[1]+coeff$zero[i])
				mu = exp(coeff$count[1]+coeff$zero[i])
				r = rzinbinom(100000,mu=mu,size=theta,zprob=zprob)
				q.mod[i,] = quantile(r,prob=c(0.95,0.99))
			}
		}
		else {
			if (i==1) {
				zprob = plogis(coeff$zero[1]+coeff$zero[5])
				mu = exp(coeff$count[1]+coeff$zero[5])
				r = rzinbinom(100000,mu=mu,size=theta,zprob=zprob)
				q.mod[i,] = quantile(r,prob=c(0.95,0.99))
			}
			else {
				zprob = plogis(coeff$zero[1]+coeff$zero[i]+coeff$zero[5])
				mu = exp(coeff$count[1]+coeff$zero[i])
				r = rzinbinom(100000,mu=mu,size=theta,zprob=zprob)
				q.mod[i,] = quantile(r,prob=c(0.95,0.99))
			}
		}	
	}
	return(q.mod)
}

nbprobs(mod12.zinb,1)
nbprobs(mod1,2)




