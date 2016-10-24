itermean<-rep(0,times=50)
for(j in 1:50){
  P1a.sim<-bugs(data,inits,parameters,"HW4_1a_model.txt",
  #debug=T,
  n.chains=1,n.iter=2000,n.burnin=1000,n.thin=1,bugs.seed=floor(runif(1,1,3)))
  
  itermean[j]<-P1a.sim$mean[[1]]
}

itermean2<-rep(0,times=50)
for(j in 1:50){
  P1a.sim<-bugs(data,inits,parameters,"HW4_1a_model.txt",
  #debug=T,
  n.chains=1,n.iter=2000,n.burnin=1000,n.thin=1,bugs.seed=floor(runif(1,1,14)))
  
  itermean2[j]<-P1a.sim$mean[[1]]
}

itermean3<-rep(0,times=50)
for(j in 1:50){
  P1a.sim<-bugs(data,inits,parameters,"HW4_1a_model.txt",
  #debug=T,
  n.chains=1,n.iter=2000,n.burnin=1000,n.thin=1,bugs.seed=floor(runif(1,1,10000)))
  
  itermean3[j]<-P1a.sim$mean[[1]]
}


