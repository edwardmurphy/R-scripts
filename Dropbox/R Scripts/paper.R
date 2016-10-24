# fit weibull
library(survival)
d <- data.frame(w = rweibull(6,shape=3,scale=5000), state=1)
s <- Surv(d$w,d$state)
sr <- survreg(s~1,dist="weibull")
print(paste("shape =",1/sr$scale))
print(paste("scale =",exp(sr$coefficients[1])))


# simple MC for weibull MLE bias
iter = 10000
shape = rep(-1,times = iter)
scale = rep(-1,times = iter)
for (i in 1:iter){
	w = rweibull(5,shape=3,scale=5000)
	d <- data.frame(r = w, state=1)
	s <- Surv(d$r,d$state)
	sr <- survreg(s~1,dist="weibull")
  shape[i] = 1/sr$scale
	scale[i] = exp(sr$coefficients[1])
}

mean(shape)
mean(scale)

# repeat but larger sample 
iter = 10000
shape = rep(-1,times = iter)
scale = rep(-1,times = iter)
for (i in 1:iter){
	w = rweibull(50,shape=3,scale=5000)
	d <- data.frame(r = w, state=1)
	s <- Surv(d$r,d$state)
	sr <- survreg(s~1,dist="weibull")
  shape[i] = 1/sr$scale
	scale[i] = exp(sr$coefficients[1])
}

mean(shape)
mean(scale)

# bootstrap
w = rweibull(5,shape=3,scale=5000)
iter = 10000
shape = rep(-1,times = iter)
scale = rep(-1,times = iter)
for (i in 1:iter){
	r <- sample(w,size=length(w),replace=TRUE)	
	d <- data.frame(r = r, state=1)
	s <- Surv(d$r,d$state)
	sr <- survreg(s~1,dist="weibull")
  shape[i] = 1/sr$scale
	scale[i] = exp(sr$coefficients[1])
}

mean(shape)  ### biased!
mean(scale)

# parametric bootstrap
d <- data.frame(r = w, state=1)
s <- Surv(d$r,d$state)
sr <- survreg(s~1,dist="weibull")
shape = 1/sr$scale
scale= exp(sr$coefficients[1])

shape.boot = rep(-1,times = iter)
scale.boot = rep(-1,times = iter)

for (i in 1:iter){
		w = rweibull(6,shape=shape,scale=scale)
		d <- data.frame(r = w, state=1)
		s <- Surv(d$r,d$state)
		sr <- survreg(s~1,dist="weibull")
  	shape.boot[i] = 1/sr$scale
		scale.boot[i] = exp(sr$coefficients[1])
}

# median rank regression
# code from Symynck XVI-th ISC

mrank.observation <- function (j, f){
	r <- qbeta(0.5, j, f - j + 1);r
}

mrank.data <- function (data = NULL){
	n <- nrow(data)
	data$rank <- rank(data$ob, ties.method =
"first")
	data$rrank <- (n + 1 - data$rank)
	sdata <- data[order(data$rank), ]
	sdata$arank <- sdata$rank
	sdata$mrank <- mrank.observation(sdata$arank, n)
	data <- sdata[!is.na(sdata$mrank), ];data
	}

d <- data.frame(ob=c(149971, 70808, 133518,
145658, 175701, 50960, 126606, 82329), state=1)
d <- mrank.data(d)
# d <- mrank.data(data.frame(ob = w,state=1))

maptowb <- function (y){log(log(1/(1 - y)))}
fwb <- lm(log(d$ob) ~ maptowb(d$mrank),d)
print(summary(fwb))
print(paste(
"beta =",1/coefficients(fwb)[2]))
print(paste(
"eta =", exp(coefficients(fwb)[1])))


