assoc_Zi.Zj <-
function(x){
len<-ncol(x)
res <- matrix(0, ncol=ncol(x), nrow=ncol(x))
  for(i in 1:len){
    for(j in 1:len){
      if(i != j) {
        sumpair <- sum(x[,i]=="C"&x[,j]=="C" | x[,i]=="C"&x[,j]=="T" |  x[,i]=="T"&x[,j]=="C" | x[,i]=="T"&x[,j]=="T")
        res[i,j] <- sumpair/nrow(x)
      }
    }
  }
p <- res

cal<-matrix(0,nrow=3,ncol=ncol(x))
  for(j in 1:len){
    nA<-sum(x[,j]=="A")
    nT<-sum(x[,j]=="T")
    nG<-sum(x[,j]=="G")
    nC<-sum(x[,j]=="C")
    cal[1,j]<-qnorm((nA+nG)/nrow(x))
    cal[2,j]<-qnorm(nT/(nC+nT))
    cal[3,j]<-qnorm(nG/(nA+nG))
     
  }
vry <- cal
#################################################################
l<-matrix(0,nrow=len,ncol=len)
for(i in 1:len) {
                      for(j in 1:len){
                      set.seed(i*j)
			if(i!=j){
                         fun<-function(s){pmvnorm(mean=c(0,0),
                        lower=c(vry[1,i],vry[1,j]),upper=c(Inf,Inf),
                        sigma=matrix(c(1,s,s,1),ncol=2))-p[i,j]}
                     
                       l[i,j]<-try(uniroot(fun,lower=-1,upper=1,extendInt="yes")$root)
                       options(show.error.messages=FALSE)
                                  
                  }
         }
       }
   diag(l) <- 1
   return(round(l, 3))
}
