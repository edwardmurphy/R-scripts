###################################################
### chunk number 1: 
###################################################
#line 16 "HW7.Rnw"
# read data
law = read.table("C:\\Users\\Owner\\Dropbox\\MATH352\\law_data.txt",header=F)
L = law[,1]
G = law[,2]

# correlation
Lbar = mean(L)
Gbar = mean(G)
r.pop = sum((L-Lbar)*(G-Gbar))/sqrt( sum( (L-Lbar)^2)*sum((G-Gbar)^2))


###################################################
### chunk number 2: 
###################################################
#line 29 "HW7.Rnw"
r.pop


###################################################
### chunk number 3: 
###################################################
#line 33 "HW7.Rnw"
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


###################################################
### chunk number 4: 
###################################################
#line 55 "HW7.Rnw"
data = rnorm(20,mean=5,sd=2)
n = length(data)
alpha = 0.05
low.t.ci = mean(data) - qt(1-alpha,n-1)*sd(data)/sqrt(n)
upper.t.ci = mean(data) + qt(1-alpha,n-1)*sd(data)/sqrt(n)

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


###################################################
### chunk number 5: t_pivot
###################################################
#line 87 "HW7.Rnw"
png("t_pivot.png")
par(mfrow=c(1,1))
hist(t.star,breaks=100,freq=F,main="",xlab="Pivot Value")
x = seq(-8,8,by=0.001)
y = dt(x,df=19)
lines(x,y)
null<-dev.off()


###################################################
### chunk number 6: 
###################################################
#line 109 "HW7.Rnw"
n = 20

#### monte carlo t pivot for exp(mean=5)
iter = 100000 # number of monte carlo iterations
t = vector("numeric",iter) # array to hold pivot/iteration

for (i in 1:iter){
	data = rexp(n,rate=1/5)
	t[i] = sqrt(n)*(mean(data)-5)/sd(data)
}


###################################################
### chunk number 7: exp_pivot
###################################################
#line 122 "HW7.Rnw"
png("exp_pivot.png")
par(mfrow=c(1,1))
hist(t,breaks=100,xlim=c(-10,4),main="",xlab="Pivot Value")
null<-dev.off()


###################################################
### chunk number 8: 
###################################################
#line 143 "HW7.Rnw"
data = rexp(n,rate=1/5)  # generate fixed data for bootstrap

B = 100000  # number of bootstrap iterations
t_star = vector("numeric",B)  # array to hold pivot/iteration

for (i in 1:B){
	data_star = sample(data,replace=T)
	t_star[i] = sqrt(n)*(mean(data_star)-mean(data))/sd(data_star)
}


###################################################
### chunk number 9: exp_boot_pivot
###################################################
#line 155 "HW7.Rnw"
png("exp_boot_pivot.png")
hist(t_star,breaks=100,xlim=c(-10,4),main="",xlab="Pivot Value",freq=F)
x = seq(-10,4,by=0.001)
y = dnorm(x)
lines(x,y)
null<-dev.off()


