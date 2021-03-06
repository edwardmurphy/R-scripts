#####################################################
# Function to get acceptance number for AQL sampling
# Inputs: p = probability 
#					batch = total batch size
#					k = sample size
#					aql = aql level
#	Outputs: Ac = acceptance number

# hypergeometric
getAccept <- function(p,batch,k,aql){
	m = floor(aql*batch/100)
	n = batch - m
	Ac = qhyper(p,m,n,k)
	return(Ac)
}



#####################################################

batch_min <- c(2,9,16,26,51,91,151,281,501,1201,3201,10001,35001,150001)
batch_max <- c(batch_min[2:14]-1,500000)
ss_min <- c(batch_min[1:7],200,200,200,200,800,800,800) 
ss_max <- c(batch_max[1:6],200,200,200,200,200,800,800,800) 

accept_min <- rep(-1,length(batch_min))
#reject_min <- rep(-1,length(batch_min))
accept_max <- rep(-1,length(batch_max))
#reject_max <- rep(-1,length(batch_max))

aql_accept_min <- matrix(rep(accept_min,3),ncol=3,byrow=F)
aql_accept_max <- matrix(rep(accept_max,3),ncol=3,byrow=F)

p <- 0.95
aql <- c(0.10,1.0,4.0)

for (j in 1:length(aql)){
	for (i in 1:length(batch_min)){
			accept_min[i] <- getAccept(p,batch_min[i],ss_min[i],aql[j])
			accept_max[i] <- getAccept(p,batch_max[i],ss_max[i],aql[j])

			#reject_min <- accept_min + 1
			#reject_max <- accept_max + 1

	}

	aql_accept_min[,j] <- cbind(accept_min)
	aql_accept_max[,j] <- cbind(accept_max)

}










# write function to get p for given sample size,accept/reject
