## load workspace

load("  ") # enter location, or use GUI (File > Load Workspace...)

## used w = 1
## sigma = I
## may not be optimal parameters, timing study only

## total time (in hrs.) for 2500 obs. in 10-space
time[3]/3600

# total time is segmented into:
#	height.time = time to evaluate the current parameter values
#		using f1x, then draw uniform random for height
#		step 1 of slice sampler
#	sample.space.time = time to construct the hypercube s.t.
#		f1x at all vertices is < height
#	M.time = time to run the optim function in R to find the
#		max of the accept-reject ratio subject to boundary
#		identified by sample.space; function is called at each
#		vertex
#	AR.time = time to run the accept-reject sampler using m.v. t-dist 
#		and also time to check conditions of slice sampler (that 
#		f1x(candidate) > height and that candidate is within 
#		hypercube 


# total time for each iteration
iter.time<-height.time+sample.space.time+M.time+AR.time

# plot time for each iteration
par(mfrow=c(2,1))
hist(iter.time,breaks=100)
hist(iter.time,breaks=100,ylim=c(0,50))

# one (majority) set takes 0-15 seconds, other set takes 20-50 seconds

# determine proportions

sum(height.time)/time[3]*100

sum(sample.space.time)/time[3]*100

sum(M.time)/time[3]*100

sum(AR.time)/time[3]*100

# M.time is 90% of total time
# sample.space.time is 9% of total time

# sample stats

apply(X,1,mean)
apply(X,1,var)

# expect variance to change in different spaces? was ~3-4 in 3 space

# plot samples 
par(mfrow=c(2,5))
for (i in 1:nrow(X)){
	hist(X[i,],breaks=50, main="")
}



