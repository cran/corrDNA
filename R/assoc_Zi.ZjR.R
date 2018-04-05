assoc_Zi.ZjR <-
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
        sumpair <- sum((x[,i]=="A" & x[,j]=="A") | (x[,i]=="G" & x[,j]=="A"))
        div<-sum((x[,j]=="A") | (x[,j]=="G"))
        res[i,j] <- sumpair/div
      }
    }
  }
 pIJr <- res  
######################################################
r <- rZiZj
###########################################################
 
res<-matrix(0,nrow=len,ncol=len)
    for(i in 1:len) {
    for(j in (1:len)){
        if(i!=j){
	set.seed(i*j)

        fun<-function(s){pmvnorm(mean=c(0,0,0),lower=c(-Inf,vry[3,j],-Inf),upper=c(vry[1,i],Inf,vry[1,j]),
               sigma=matrix(c(1,s,r[i,j],s,1,0,r[i,j],0,1),ncol=3))- (pnorm(vry[1,j])*pIJr[i,j])}
        
        res[i,j]<-try(uniroot(fun,lower=-sqrt(1-(r[i,j]^2)),upper=sqrt(1-(r[i,j]^2)),extendInt="yes")$root)
		 options(show.error.messages=FALSE)
	  }  
    }
    
  }
  diag(res)<-0
  return(round(res, 3))

}
