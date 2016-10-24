x1<-runif(10000,-1,1)
x2<-runif(10000,-1,1)

fx<-rep(0,length(x1))
for (i in 1:nrow(supp.rand)){

	fx[i]<- sqrt(1-x1[i]^2-x2[i]^2) 
}

top<-scatterplot3d(x1,x2,fx,box=F,angle=24, pch=20,zlim=c(-1,1))
top$points3d(x1,x2,-fx)