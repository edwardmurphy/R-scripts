############### 
# Problem 1
###############

# read data
law = read.table("C:\\Users\\emurphy\\Dropbox\\MATH352\\law_data.txt",header=F)
law = read.table("C:\\Users\\Owner\\Dropbox\\MATH352\\law_data.txt",header=F)
L = law[,1]
G = law[,2]

# correlation
Lbar = mean(L)
Gbar = mean(G)
r.pop = sum((L-Lbar)*(G-Gbar))/sqrt( sum( (L-Lbar)^2)*sum((G-Gbar)^2))

# get sample data ("observed")
indices = c(4,6,13,15,31,35,36,45,47,50,52,53,70,79,82)
smp.data = law[indices,]

# bootstrap
B = 200
rstar = vector("numeric",length=B)

for(i in 1:B){
	data.star = sample(smp.data,length(smp.data),replace=T)
	Lbar.star = mean(data.star[,1])
	Gbar.star = mean(data.star[,2])
	rstar[i] = sum((data.star[,1]-Lbar.star)*(data.star[,2]-Gbar.star))/
		sqrt( sum( (data.star[,1]-Lbar.star)^2)*sum((data.star[,2]-Gbar.star)^2))
}

sd(rstar)

cat("The bootstrap estimate of SE(r) is",round(sd(rstar),3),"\n")

############### 
# Problem 2
###############

# create data
data = rnorm(20,mean=5,sd=2)
n = length(data)
alpha = 0.05
low.t.ci = mean(data) - qt(1-alpha,n-1)*sd(data)/sqrt(n)
upper.t.ci = mean(data) + qt(1-alpha,n-1)*sd(data)/sqrt(n)

# bootstrap
B = 10000
t.star = vector("numeric",B)

for (i in 1:B){
	data.star = sample(data,n,replace=TRUE)
	t.star[i] = sqrt(n)*(mean(data.star)-mean(data))/sd(data.star)
}

a=quantile(t.star,prob=0.05)
b=quantile(t.star,prob=0.95)

low.boot.ci = mean(data) - b*sd(data)/sqrt(n)
upper.boot.ci = mean(data) - a*sd(data)/sqrt(n)


# plot results and compare to true t
par(mfrow=c(1,1))
hist(t.star,breaks=100,freq=F)
x = seq(-8,8,by=0.001)
y = dt(x,df=19)
lines(x,y)

# get quantiles, symmetric?

quantile(t.star,probs=c(alpha,(1-alpha)))


############### 
# Problem 3
###############
n = 20

#### monte carlo t pivot for exp(mean=5)
iter = 100000 # number of monte carlo iterations
t = vector("numeric",iter) # array to hold pivot/iteration

for (i in 1:iter){
	data = rexp(n,rate=1/5)
	t[i] = sqrt(n)*(mean(data)-5)/sd(data)
}

#### bootstrap t pivot for given random data from exp(mean=5)
data = rexp(n,rate=1/5)  # generate fixed data for bootstrap

B = 100000  # number of bootstrap iterations
t_star = vector("numeric",B)  # array to hold pivot/iteration

for (i in 1:B){
	data_star = sample(data,replace=T)
	t_star[i] = sqrt(n)*(mean(data_star)-mean(data))/sd(data_star)
}

par(mfrow=c(2,1))
hist(t,breaks=100,xlim=c(-10,4))
hist(t_star,breaks=100,xlim=c(-10,4))

# compare 

cat("The t-based 90% CI is",quantile(t,prob=c(0.025,0.975)),"\n")
cat("The bootstrap 90% CI is",quantile(t_star,prob=c(0.025,0.975)),"\n")