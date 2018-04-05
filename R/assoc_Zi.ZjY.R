assoc_Zi.ZjY <-
function(x, rZiZj){
  len<-ncol(x)
  
###########################################  
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
#################################################################################  
 res <- matrix(0, ncol=ncol(x), nrow=ncol(x))
  for(i in 1:len){
    for(j in 1:len){
      if(i != j) {
        sumpair <- sum(x[,i]=="C"&x[,j]=="C" | x[,i]=="T"&x[,j]=="C")
        div<-sum(x[,j]=="C",x[,j]=="T")
        res[i,j] <- sumpair/div
      }
    }
  }
 pIJy <- res  
######################################################
r <- rZiZj
###########################################################
 
l<-matrix(0,nrow=len,ncol=len)
  for(i in 1:len) {
    for(j in 1:len){
      if(i!=j){
      set.seed(i*j)
        fun<-function(s){pmvnorm(mean=c(0,0,0),
                          lower=c(vry[1,i],vry[2,j],vry[1,j]),upper=c(Inf,Inf,Inf),
                          sigma=matrix(c(1,s,r[i,j],s,1,0,r[i,j],0,1),ncol=3))-((1-pnorm(vry[1,j]))*pIJy[i,j])}
        
        l[i,j]<-try(uniroot(fun,lower=-sqrt(1-(r[i,j]^2)),upper=sqrt(1-(r[i,j]^2)),extendInt="yes")$root)
		 options(show.error.messages=FALSE)
        
       }
    }
    
  }
  diag(l)<-0
  return(round(l, 3))

}
