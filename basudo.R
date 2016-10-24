#Attempt at a simulated annealing sudoku solver

#Sudoku sheet: 0 means empty

#Reference value (hard sudoku)
#s=matrix(0,ncol=9,nrow=9)
#s[1,8:9]=c(3,9)
#s[2,c(6,9)]=c(1,5)
#s[3,c(3,5,7)]=c(3,5,8)
#s[4,c(3,5,9)]=c(8,9,6)
#s[5,c(2,6)]=c(7,2)
#s[6,c(1,4)]=c(1,4)
#s[7,c(3,5,8)]=c(9,8,5)
#s[8,c(2,7)]=c(2,6)
#s[9,c(1,4)]=c(4,7)

#Too easy one 
#s=matrix(0,ncol=9,nrow=9)
#s[1,c(6,8)]=c(6,4)
##s[2,c(1:3,8)]=c(2,7,9,5)
#s[3,c(2,4,9)]=c(5,8,2)
#s[4,3:4]=c(2,6)
#s[6,c(3,5,7:9)]=c(1,9,6,7,3)
#s[7,c(1,3:4,7)]=c(8,5,2,4)
#s[8,c(1,8:9)]=c(3,8,5)
#s[9,c(1,7,9)]=c(6,9,1)

#Le Monde 241008
#s=matrix(0,ncol=9,nrow=9)
#s[1,c(1,5,6,8)]=c(4,9,7,8)
#s[2,2]=2
#s[3,c(3,8)]=c(6,4)
#s[4,9]=9
#s[5,c(2,4,7)]=c(1,4,6)
#s[6,c(8,3,6,8,9)]=c(8,3,2,5,7)
#s[7,1:3]=c(5,9,7)
#s[8,5]=1
#s[9,5:7]=c(2,6,7)

#Le Monde 210210
s=matrix(0,ncol=9,nrow=9)
s[1,c(1,6,7)]=c(8,1,2)
s[2,c(2:3)]=c(7,5)
s[3,c(5,8,9)]=c(5,6,4)
s[4,c(3,9)]=c(7,6)
s[5,c(1,4)]=c(9,7)
s[6,c(1,2,6,8,9)]=c(5,2,9,4,7)
s[7,c(1:3)]=c(2,3,1)
s[8,c(3,5,7,9)]=c(6,2,1,9)

# target function to minimise
# that counts lines+rows+blocks
# missing the no duplicate constraint
target=function(s){

  tar=sum(apply(s,1,duplicated)+apply(s,2,duplicated))

   for (r in 1:9){

      bloa=(1:3)+3*(r-1)%%3
      blob=(1:3)+3*trunc((r-1)/3)
      tar=tar+sum(duplicated(as.vector(s[bloa,blob])))

     }

  return(tar)
  } 

#local score
scor=function(i,s){

 a=((i-1)%%9)+1
 b=trunc((i-1)/9)
 boxa=3*trunc((a-1)/3)+1
 boxb=3*trunc(b/3)+1

 return(sum(s[i]==s[9*b+(1:9)])+
      sum(s[i]==s[a+9*(0:8)])+
      sum(s[i]==s[boxa:(boxa+2),boxb:(boxb+2)])-3)
 }

#Deterministic prunning
pool=array(TRUE,dim=c(9,9,9))

for (i in 1:9)
for (j in 1:9){

  if (s[i,j]>0) pool[i,j,-s[i,j]]=FALSE
  }

for (t in 1:100){ # random order for visit of all sites

  for (i in sample(1:81)){

    if (s[i]==0){
 
       a=((i-1)%%9)+1
       b=trunc((i-1)/9)+1
       boxa=3*trunc((a-1)/3)+1
       boxa=boxa:(boxa+2)
       boxb=3*trunc((b-1)/3)+1
       boxb=boxb:(boxb+2)

       for (u in (1:9)[pool[a,b,]]){#eliminates impossible values

         pool[a,b,u]=(sum(u==s[a,])+sum(u==s[,b])
		+sum(u==s[boxa,boxb]))==0
         }

       if (sum(pool[a,b,])==1){ # only one possible case, solution found!

         s[i]=(1:9)[pool[a,b,]]
         print(s)
         }

       if (sum(pool[a,b,])==0){ # solution does not exist, exit!

	  print("wrong sudoku")
          break()
          }
   
       }
    }
  }

#Simulated annealling sequence

#Parameters
#Uneducated guesses
Niter=10^4
Nsteps=10^2
lmax=10^5
temps=exp(sqrt(seq(1,log(lmax)^2,le=Niter+1)))
lcur=temps[1]

cur=s
for (r in (1:81)[s==0])
     cur[r]=sample((1:9)[pool[r+ 81*(0:8)]],1)
tarcur=target(cur)

bess=cur
bestar=tarcur

plot(c(0,0),col="white",xlim=c(1,Niter),ylim=c(0,tarcur),xlab="iterations",ylab="penalty")

for (t in 1:Niter){

  nchange=0
  aveva=0*s
  naver=0

  for (d in 1:Nsteps){ 

    prop=cur

    # random moves on the sites
    i=sample((1:81)[as.vector(s)==0],sample(1:sum(s==0),1,pro=1/(1:sum(s==0))))
    for (r in 1:length(i))
      prop[i[r]]=sample((1:9)[pool[i[r]+81*(0:8)]],1)

    if (log(runif(1))/lcur<tarcur-target(prop)){

        nchange=nchange+(tarcur>target(prop))
	cur=prop
points(t,tarcur,col="forestgreen",cex=.3,pch=19)
        tarcur=target(cur)
        }

    if (tarcur==0)
	break()

    aveva=aveva+cur
    naver=naver+1
    if (tarcur<bestar){

	bestar=tarcur
	bess=cur
	}

    lcur=sample(c(1,10^(-4)),1,pro=c(1.5*(log(t+1))+1,1))*temps[t+1]

    if (tarcur==0)
	break()

    #Another basic but efficient go
    for (r in sample((1:81)[s==0])){

      prop=cur

      if (scor(r,cur)>0){

        prop[r]=sample((1:9)[pool[r+81*(0:8)]],1)

        if (log(runif(1))/lcur<tarcur-target(prop)){

          nchange=nchange+(tarcur>target(prop))
          cur=prop
          points(t,tarcur,col="gold",cex=.3,pch=19)
          tarcur=target(cur)
          }
        }

      }
} # end of the t loop
