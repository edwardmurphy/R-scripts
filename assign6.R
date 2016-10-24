### Stroke data
s1 = 119
s2 = 98
n1 = 11037
n2 = 11034

array1 = c(rep(1,s1),rep(0,n1-s1))
array2 = c(rep(1,s2),rep(0,n2-s2))

B = 1000

ratio.star = vector("numeric",B)

for (i in 1:B){
	s1.star = sum(sample(array1,n1,replace=TRUE))
	s2.star = sum(sample(array2,n2,replace=TRUE))
	ratio.star[i] = (s1.star/n1)/(s2.star/n2)
}

hist(ratio.star,breaks=50)
quantile(ratio.star,probs=c(0.025,0.975))

### Normal example
data = rnorm(20,mean=5,sd=2)
n = length(data)

#normal CI
norm.low.ci <- mean(data) - qnorm(0.95)*sd(data)/sqrt(n)
norm.up.ci <- mean(data) + qnorm(0.95)*sd(data)/sqrt(n)

#t CI
t.low.ci <- mean(data) - qt(0.95,n-1)*sd(data)/sqrt(n)
t.up.ci <- mean(data) + qt(0.95,n-1)*sd(data)/sqrt(n)

#bootstrap
B = 10000
mean.star = vector("numeric",B)

for (i in 1:B){
	mean.star[i] = mean(sample(data,n,replace=TRUE))
}

boot.ci <- quantile(mean.star,probs=c(0.05,0.95))


