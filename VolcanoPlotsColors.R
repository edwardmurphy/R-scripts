# --------------------------------------------------
# Math 365A Computational Methods Molecular Biology
# Prof. Claudia Rangel-Escareño
# 
# Colored Volcanoplots
# --------------------------------------------------

# First look at the results in B&W to determine the axes limits 
# to avoid misplacing the labels every time I overwrite colors
# If nC is the number of contrasts of interest (suppose 3)
 
par(mfrow=c(1,nC))
for (i in 1:nC)
{
volcanoplot(fit2, coef=i)
}

# These are vectors of as many entries as nC update the values
# accordingly 

y0=c(-5,-6,-7)
y1=c(-2,3,5)
x0=c(-3,-3,-4)
x1=c(2,4,3)

# Next, I set the values of the B statistic (B) and the log fold-change
# which I call M here. There is also another variable "contrasts"
# This is the vector of labels for each contrast "this vs that"

M=1
B=0

for (i in 1:nC)
 {
    volcanoplot(fit2, col="blue", ylim=c(y0[i],y1[i]), xlim=c(x0[i],x1[i]), coef=i, main=contrasts[i], cex.lab=1.3)
    par(new=T)
    abline(v=-M, col="brown")
    par(new=T)
    abline(v=M, col="brown")
    par(new=T)
    abline(h=B, col="black")
    par(new=T)
    ind1 = abs(fit2$coef[,i])>M
    ind2 = fit2$lods[,i] >B
    ind3 = (fit2$coef[,i]>M & fit2$lods[,i]>B)
    ind4 = (fit2$coef[,i]< -M & fit2$lods[,i]>B)
    x = as.matrix(fit2$coef[ind1,i])
    y = as.matrix(fit2$lods[ind1,i])
    plot(x, y, col="magenta",ylim=c(y0[i],y1[i]), xlim=c(x0[i],x1[i]),main="", pch = "*",cex.lab=1.2)
    x = as.matrix(fit2$coef[ind2,i])
    y = as.matrix(fit2$lods[ind2,i])
    par(new=T)
    plot(x, y, col="orange",  ylim=c(y0[i],y1[i]), xlim=c(x0[i],x1[i]), main="", pch = 19,xlab="", ylab="",cex.lab=1.2)
    x = as.matrix(fit2$coef[ind3,i])
    y = as.matrix(fit2$lods[ind3,i])
    par(new=T)
    plot(x, y, col="red",  ylim=c(y0[i],y1[i]), xlim=c(x0[i],x1[i]),main="", pch = 19,xlab="", ylab="",cex.lab=1.2)
    x = as.matrix(fit2$coef[ind4,i])
    y = as.matrix(fit2$lods[ind4,i])
    par(new=T)
    plot(x, y, col="darkgreen",  ylim=c(y0[i],y1[i]), xlim=c(x0[i],x1[i]),main="", pch = 19,xlab="", ylab="",cex.lab=1.2)
 }
