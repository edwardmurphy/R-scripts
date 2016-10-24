setwd()

library(affy)

mydata=ReadAffy()

hist(mydata)

sampleNames(mydata)
ma.labels = sampleNames(mydata)
ma.lab = substr(ma.labels, 1, nchar(ma.labels)-4)
sampleNames(mydata) = ma.lab
sampleNames(mydata)


image(mydata[,1])

boxplot(mydata)

mycolors = rep(c("bisque","green","salmon"), each=3)
mycolors = rep(c(bisque'',''green'',''salmon''), 3)
boxplot(mydata, col=mycolors, main=“Raw Data”)

pms = pm(mydata[,1])
mms = mm(mydata[,1])
pm.mm = mean(pms>mms)
pm.mm = 100*pm.mm

i=1
j=2
array.i = pm(mydata[,i])
array.j = pm(mydata[,j])

M = log2(array.i) – log2(array.j)
A = 0.5*(log2(array.i) + log2(array.j))
plot(A,M)
par(new=T)
abline(h=0)

