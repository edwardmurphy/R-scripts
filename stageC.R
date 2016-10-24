## There are two tree packages for a tree analysis: tree (for classification
## and regression trees) and rpart (Recursive Partitioning). The rpart
## package is newer and also has a nice accompanying documentation - An
## Introduction to Recursive Partitioning Using the RPART Routines,
## by Atkinson and Therneau (available at the course web site). We will
## use rpart in this class.

library(rpart)

############################################################
# Read in stagec.data 
# 
#   pgtime = time to progression in years
#   pgstat = status at last follow-up
#            1=progressed, 0=censored
#   age = age at diagnosis
#   eet = early endocrine therapy: 1=no 2=yes
#   g2 = % of cells in g2 phase, from flow cytometry
#   grade = tumor grade 1,2,3,4
#   gleason = Gleason score (competing grading system, 3-10)
#   ploidy = diploid/tetraploid/aneuploid DNA pattern
############################################################

stagec <-
read.table("http://www-rohan.sdsu.edu/~jjfan/sta702/stagec.txt",
col.names=c("pgtime", "pgstat", "age", "eet", "g2", "grade",
"gleason", "ploidy")) 

#####################################################
# Specify that the variable "pgstat" and "ploidy" are
# categorical and name their categories: 
#####################################################

stagec$pgstat <- factor(stagec$pgstat, levels=0:1, labels=c("No", "Prog"))
stagec$ploidy <- factor(stagec$ploidy, levels=1:3,
			labels=c("diploid", "tetraploid", "aneuploid"))

######################################################
# fit classification tree to data in stagec data frame
# using 10-fold cross-validation; method="class"
######################################################

my.control <- rpart.control(xval=10, cp=0) 
cfit1 <- rpart(pgstat ~ age+eet+g2+grade+gleason+ploidy,
	       data=stagec, method="class", control=my.control)

########################################
# plot the tree that corresponds to cp=0
# (the largest, unpruned tree):
########################################

plot(cfit1, margin=0.1)
text(cfit1, use.n=T)

########################################################
# look at the sequence of unpruned trees to decide which
# value of cp gives the optimal tree:
########################################################

printcp(cfit1)
plotcp(cfit1)

#################################################################
# plot the cross-validation estimates of error and training error
# against the tree complexity with error bars:
#################################################################

numsplits <- cfit1$cptable[, 2]
trainerror <- cfit1$cptable[, 3]
xerror <- cfit1$cptable[, 4]
xstd <- cfit1$cptable[, 5]

plot(numsplits, trainerror, ylim=c(0.5, 1.2), type="l")
lines(numsplits, xerror, lty=2)
lines(numsplits, xerror-xstd, lty=2)
lines(numsplits, xerror+xstd, lty=2)
title("Cross-validation Error Estimates and Training Error")
legend(.02, .55, c("trainerror", "xerror"), lty=c(1,2))

#################################
# get the optimal tree by pruning
#################################

cfit1pruned <- prune(cfit1, cp=.03)

print(cfit1pruned)
plot(cfit1pruned, uniform=T)
text(cfit1pruned, use.n=T)
summary(cfit1pruned)
summary(cfit1, cp=0.03)  # same as above

#####################################
# try different priors on the classes
#####################################

cfit2 <- rpart(pgstat ~ age+eet+g2+grade+gleason+ploidy,
	       stagec, parm=list(prior=c(.5,.5)))

printcp(cfit2)
print(cfit2)
print(cfit2,cp=.05)
summary(cfit2)
summary(cfit2,cp=.05)

############################################################
# the following code shows various ways to plot your results
############################################################

plot(cfit1pruned,uniform=T)
text(cfit1pruned,use.n=T,all=T)

plot(cfit1pruned,branch=.5,compress=T,uniform=T)
text(cfit1pruned,digits=2,use.n=T)

plot(cfit1pruned,compress=T)
text(cfit1pruned,digits=2)

#################################################################
## post.rpart creates a postscript file with a nicer looking tree
#################################################################

post(cfit1pruned,title="Stage C prostate Cancer Data",file="",use.n=T)

