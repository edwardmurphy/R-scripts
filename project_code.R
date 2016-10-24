### add mixture model to randomly draw from module which fails?
### why? maybe modules have variable repair times and costs
### 


library(lattice)
library(xtable)

########################################################
##### Play around with Weibull distribution
##### Note that shape < 1 : DFR
##### 					shape = 1 : exponential
#####						shape > 1 : IFR
#####						shape = 2 : Rayleigh (linear hazard)

# parameters of weibull
a <- 3			#shape
b <- 1500		#scale

# moments
EX <- b * gamma(1 + 1/a)
VarX <-  b^2 * (gamma(1 + 2/a) - (gamma(1 + 1/a))^2)

# grid
x <- seq(0,EX+6*sqrt(VarX),by=0.1)

# distributions
dens <- dweibull(x,shape=a,scale=b)
cdf <- pweibull(x,shape=a,scale=b)
reli <- 1-cdf
haz <- dens/reli

plot(x,dens,type="l")
plot(x,reli,type="l")
plot(x,haz,type="l")

########################################################

########################################################
##### Play around with lognormal distribution
mu <- 3
sigma <- 1

# parameters
logNormMoments <- function(mu,sigma){
	# moments
	EX <- exp(mu + 1/2*sigma^2)
	VarX = exp(2*mu + sigma^2)*(exp(sigma^2) - 1)
	return(c(mu,sigma,EX,VarX,sqrt(VarX)))
}

ln.mom <- data.frame(c(0,0,0,0,0))
rownames(ln.mom) <- c("mu","sigma","EX","VarX","SigmaX")

ln.mom[,1] <- logNormMoments(0,1)
ln.mom[,2] <- logNormMoments(-1,1)
ln.mom[,3] <- logNormMoments(0.5,1)
ln.mom[,4] <- logNormMoments(-0.5,1)
ln.mom[,5] <- logNormMoments(0.5,0.5)
ln.mom[,6] <- logNormMoments(-0.5,0.5)
ln.mom[,7] <- logNormMoments(1,1)
ln.mom[,8] <- logNormMoments(10,.1)

### sigma must be small (> 0)- variance explodes even at sigma = 2
### (0,0) --> EX = 1 w.p. 1
### (0,1) --> EX = 1.64, SigX = 2.16
### use mu < 0 to get EX < 1

# grid
x <- seq(0,EX+6*sqrt(VarX),by=0.1)

# distributions
dens <- dlnorm(x,mu,sigma)
cdf <- plnorm(x,mu,sigma)
reli <- 1-cdf
haz <- dens/reli

plot(x,dens,type="l")
plot(x,reli,type="l")
plot(x,haz,type="l")

########################################################


########################################################
##### Simple Monte Carlo Example
running <- 0
down <- 0
availability <- 0
iter <- 10000

for (i in 1:iter){
	while (running < 1500){
		u <- rweibull(1, shape = 3, scale = 100) #mean = scale * gamma(1 + 1/shape)
		running <- running + u
		d <- rlnorm(1, 1, 0.5)
		down <- down + d
	}
	availability[i] <- running/(running+down)
	running <- 0
	down <- 0
}

mean(availability)
########################################################

########################################################
##### Simple Monte Carlo Example with Dependent Repair
running <- 0
down <- 0
availability <- 0
iter <- 10000

for (i in 1:iter){
	while (running < 1500){
		u <- rweibull(1, shape = 3, scale = 100) #mean = scale * gamma(1 + 1/shape)
		running <- running + u
		if (u < 70) c <- 2 * runif(1)	else c <-1
		d <- rlnorm(1, c*1, c*0.5)
		down <- down + d
	}
	availability[i] <- running/(running+down)
	running <- 0
	down <- 0
}

mean(availability)
########################################################


########################################################
##### Compare algorithms (quantile function vs direct)
##### 1. Draw uniform random, determine time
##### 2. Draw random time

n<-50000
# 1
u1<-runif(n)
w1<-qweibull(u,shape=2,scale=10000)

# 2
w2<-rweibull(n,shape=2,scale=10000)

# compare
par(mfrow=c(2,1))
hist(w1,breaks=100)
hist(w2,breaks=100)
########################################################


########################################################
##### Simulate from different failure distributions
##### log availability

current.up.time <- 0
total.up.hours <- 20000

# weibull shape and scale parameters for failure distribution
fail.params <- matrix(c(2,1000,10,5000,20,10000),nrow=2,byrow=F)

# parameters of repair distribution
repair.params <- matrix(c(0,1,0,0.5,2,0.25),nrow=2,byrow=F)

# initialize list of failure times for each component
fail.times <- matrix(rep(0,times=3),nrow=1)

fail.times[1,1] <- rweibull(1,shape=fail.params[1,1],scale=fail.params[2,1])
fail.times[1,2] <- rweibull(1,shape=fail.params[1,2],scale=fail.params[2,2])
fail.times[1,3] <- rweibull(1,shape=fail.params[1,3],scale=fail.params[2,3])

# initialize list of component failure times
# will collect failure times for each component
comp.fail.times <- list(0,0,0)

# initialize list of component repair times
# will collect repair times for each component
comp.repair.times <- list(0,0,0)

# event-driven simulation
while (current.up.time < total.up.hours){    	#run sim for certain up hours
	new.fail.time <- min(fail.times)						#new fail time is min time in row
	j <- which.min(fail.times)									#comp (col) which failed
	
	# save fail time for component
	comp.fail.times[[j]] <- rbind(comp.fail.times[[j]],new.fail.time)

	# sample down time for failed component
	comp.repair.times[[j]] <- rbind(comp.repair.times[[j]],rlnorm(1,repair.params[1,j],repair.params[2,j]))
	
	# resample for failed component
	# note that this leaves the times for components which have not failed
	fail.times[1,j] <- new.fail.time + rweibull(1,shape=fail.params[1,j], scale=fail.params[2,j])
		
	current.up.time <- new.fail.time						#update current up time
}

# compute time between failures
comp.fail.rates <- list(0,0,0)
for (j in 1:length(comp.fail.times)){
	for (k in 2:length(comp.fail.times[[j]])){
			comp.fail.rates[[j]][k-1]<-comp.fail.times[[j]][k]-comp.fail.times[[j]][k-1]
	
	}
}
	
# compute availability for each component
comp.avail <- list(0,0,0)
comp.sum.repair.times <- lapply(comp.repair.times,sum)
for (i in 1:length(comp.repair.times)){
	comp.avail[[i]] <- current.up.time/(current.up.time+comp.sum.repair.times[[i]])
}

# compute system availability
sys.sum.repair.times <- sum(unlist(comp.sum.repair.times))
sys.avail <- current.up.time/(current.up.time+sys.sum.repair.times)

# plot failure times
stripchart(comp.fail.times[[1]][-1,1],pch=1)
stripchart(comp.fail.times[[2]][-1,1],pch=2,add=T)
stripchart(comp.fail.times[[3]][-1,1],pch=3,add=T)


########################################################
##### Table with failure/repair parameters for Agilent

# Pump seals
# Pistons
# AIV
# Outlet valve
# Filter frit
# Degasser

# Needle 
# Needle seat
# Needle arm
# Chiller
# Lamp
# Flow cell
########################################################


########################################################
##### Availability/cost analysis WITHOUT PM
# initialize
iter <- 1000

comps <- c("pumpseal","piston1","piston2","inletvalve","outletvalve",
						"filterfrit","needle","needleseat","needlearm","lamp","flowcell")

# weibull shape and scale parameters for failure distribution
fail.params <- matrix(c(3,2000,  #pump seals
												5,3000,  # piston 1
												5,3000,  # piston 2
												5,3000,  # inlet valve
												5,3000,  # outlet valve
												3,2000,  # filter frit
												5,3000,  # needle 
												4,2000,  # needle seat
												20,10000, # needle arm
												3,1500,  # lamp
												10,3000), # flow cell
												nrow=2,byrow=F)

colnames(fail.params) <- comps

# initialize list of failure times for each component
fail.times <- matrix(rep(0,times=11),nrow=1)
colnames(fail.times) <- comps

# parameters of repair distribution
repair.params <- matrix(c(0,1,  		#pump seals
												0,1,  			# piston 1
												0,1,  			# piston 2
												1,0.25,  		# inlet valve
												1,0.25,  		# outlet valve
												-0.5,0.1,  	# filter frit
												1,0.25,  		# needle 
												1,0.25,  		# needle seat
												3,1,			 	# needle arm
												1,0.25,		  # lamp
												1,0.25), # flow cell
												nrow=2,byrow=F)

colnames(repair.params) <- comps

comp.avail <- list(0,0,0,0,0,0,0,0,0,0,0)
names(comp.avail) <- comps
sys.avail <- rep(0,times=iter)

comp.cost <- list(60,150,150,200,200,15,75,45,2000,800,1500)
	names(comp.cost) <- comps
	comp.repair.cost <- list(0,0,0,0,0,0,0,0,0,0,0)

## ITERATION
for (I in 1:iter){

	current.up.time <- 0
	total.up.hours <- 2184*4*5  # sim for 5 years (each PM = 2184 hrs)

	for (i in 1:length(fail.times)){
	fail.times[1,i] <- rweibull(1,shape=fail.params[1,i],scale=fail.params[2,i])
	}

	# initialize list of component failure times
	# will collect failure times for each component
	comp.fail.times <- list(0,0,0,0,0,0,0,0,0,0,0)
	names(comp.fail.times) <- comps

	# initialize list of component repair times
	# will collect repair times for each component
	comp.repair.times <- list(0,0,0,0,0,0,0,0,0,0,0)
	names(comp.repair.times) <- comps

	# event-driven simulation
	while (current.up.time < total.up.hours){    	#run sim for certain up hours
		new.fail.time <- min(fail.times)						#new fail time is min time in row
		j <- which.min(fail.times)									#comp (col) which failed
	
		# save fail time for component
		comp.fail.times[[j]] <- rbind(comp.fail.times[[j]],new.fail.time)

		# sample down time for failed component
		comp.repair.times[[j]] <- rbind(comp.repair.times[[j]],rlnorm(1,repair.params[1,j],repair.params[2,j]))
	
		# resample for failed component
		# note that this leaves the times for components which have not failed
		fail.times[1,j] <- new.fail.time + rweibull(1,shape=fail.params[1,j], scale=fail.params[2,j])
		
		current.up.time <- new.fail.time						#update current up time
	}

	# compute time between failures
	comp.fail.rates <- list(0,0,0,0,0,0,0,0,0,0,0)
	names(comp.fail.rates) <- comps
	for (j in 1:length(comp.fail.times)){
		for (k in 2:length(comp.fail.times[[j]])){
				comp.fail.rates[[j]][k-1]<-comp.fail.times[[j]][k]-comp.fail.times[[j]][k-1]
	
		}
	}
	
	# compute availability for each component
	comp.sum.repair.times <- lapply(comp.repair.times,sum)
	for (i in 1:length(comp.repair.times)){
		comp.avail[[i]][I] <- current.up.time/(current.up.time+comp.sum.repair.times[[i]])
	}

	# compute system availability
	sys.sum.repair.times <- sum(unlist(comp.sum.repair.times))
	sys.avail[I] <- current.up.time/(current.up.time+sys.sum.repair.times)

	# cost analysis
	for (i in 1:length(comp.cost)){
		comp.repair.cost[[i]][I] <- comp.cost[[i]] * length(comp.fail.times[[i]])	
	}


}

comp.avail.mean <- lapply(comp.avail,mean)
sys.avail.mean <- mean(sys.avail)
comp.repair.cost.mean <- lapply(comp.repair.cost,mean)
sys.repair.cost.year.mean <- sum(unlist(comp.repair.cost.mean))/5


sys.repair.cost.life <-mapply(sum,comp.repair.cost[[1]],comp.repair.cost[[2]],comp.repair.cost[[3]],comp.repair.cost[[4]],comp.repair.cost[[5]],comp.repair.cost[[6]],comp.repair.cost[[7]],comp.repair.cost[[8]],comp.repair.cost[[9]],comp.repair.cost[[10]],comp.repair.cost[[11]])

sys.repair.cost.year <- sys.repair.cost.life/5
 
availTable <- cbind(unlist(comp.avail.mean),unlist(lapply(comp.avail,min)),unlist(lapply(comp.avail,max)))

sysVec <- c(sys.avail.mean,min(sys.avail),max(sys.avail))

costVec <- c(mean(sys.repair.cost.year),min(sys.repair.cost.year),max(sys.repair.cost.year))

availTable <- rbind(availTable,sysVec,costVec)

xtable(availTable,digits=5)







