a<-0
n<-1

while(n<10000){
	u<-runif(8)
	if(max(u)>0.98) a<-a+1
	n<-n+1
}

a<-0
n<-1

while(n<10000){
	u<-runif(8)
	if(min(u)>0.044) a<-a+1
	n<-n+1
}

