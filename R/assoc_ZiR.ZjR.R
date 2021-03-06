assoc_ZiR.ZjR <-
function(x,rZiZj,rZiZjR){

    len<-ncol(x)
#######################################################
r <- rZiZj
rIrJ <- t(rZiZjR) 
rIJr <- rZiZjR
########################################################
my.uniroot <- function(fun, low, high)
{	
	l <- NA
	res <- try(uniroot(fun,lower=low , upper= high))
	options(show.error.messages=FALSE)
	if(!inherits(res, "try-error")) { l <- res$root ; return(l) ; }

	ss <- seq(low, high, length=10)
	fss <- sapply(ss, fun)
	if(any(fss == 0)) { l <- ss[which(fss == 0)[1]] ; return(l) ; }
	vv <- fss[1 : (length(ss) -1)] * fss[2 : length(ss)]
	if(any(vv < 0))
	{
		ww <- which(vv < 0)
		mm <- pmin(abs(ss[ww] - low), abs(ss[ww + 1] - high))
		pp <- which.max(mm)
		low <- ss[ww[pp]]
		high <- ss[ww[pp] + 1]
		res <- try(uniroot(fun,lower=low , upper= high,extendInt="yes"))
		options(show.error.messages=FALSE)
		if(!inherits(res, "try-error")) l <- res$root
		#else print(c(low, high, ss[ww[pp]], ss[ww[pp]+1], fun(ss[pp]), fun(ss[ww[pp]]+1)))

	} else
	{
		ss1 <- seq(low, high, length=100)
		fss1 <- sapply(ss1, fun)
		l <- ss1[which.min(abs(fss1))]
	}
	l
}


##########################################################  
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
  
###########################################################
res <- matrix(0, ncol=ncol(x), nrow=ncol(x))
  for(i in 1:len){
    for(j in 1:len){
      if(i != j) {
        sumpair <- sum(x[,i]=="A"&x[,j]=="A")
	div<-sum(x[,i]=="A"&x[,j]=="A" | x[,i]=="A"&x[,j]=="G" | x[,i]=="G"&x[,j]=="A" | x[,i]=="G"&x[,j]=="G")
        res[i,j] <- sumpair/div
      }
    }
  }
  
pIrJr <- res
######################################################  
  l<-matrix(0,nrow=len,ncol=len)
  for(i in 1:len) {
    for(j in i:len){
	set.seed(i*j)
	
     
	
        fun<-function(x){pmvnorm(mean=c(0,0,0,0), lower=c(vry[3,i],vry[3,j],-Inf,-Inf), upper=c(Inf,Inf,vry[1,i],vry[1,j]),
                sigma=matrix(c(1,x,0,rIrJ[i,j],x,1,rIrJ[j,i],0,0,rIJr[i,j],1,r[i,j],rIJr[j,i],0,r[i,j],1),ncol=4))-
        ((pmvnorm(mean=c(0,0),lower=c(-Inf,-Inf),upper=c(vry[1,i],vry[1,j]),sigma=matrix(c(1,r[i,j],r[i,j],1),ncol=2)))*pIrJr[i,j])}
        
        a<-(1-(r[i,j]^2))

         b<-2 * (rIrJ[i,j]*rIJr[i,j]*r[i,j])

         c<-((rIrJ[i,j]^2)+(rIJr[i,j]^2)+(r[i,j]^2)-(rIrJ[i,j]*rIJr[i,j])^2-1)

	
         low3 <- -sqrt(1 - rIJr[i,j]^2)
	 high3 <- sqrt(1 - rIJr[i,j]^2)	
	 low2 <- -sqrt(1 - rIrJ[i,j]^2)
	 high2 <- sqrt(1 - rIrJ[i,j]^2)
	 delta <- (b^2-(4*a*c))
	 if(abs(delta) < 1e-10) delta <- 0	
             if(delta >=0){
		         		low1 <- ((-b)-sqrt(delta))/2*a
	 				high1 <- ((-b) + sqrt(delta))/2*a
					low <- max(low1, low2, low3) 
					high <- min(high1, high2, high3)					
					l[i, j] <-try(my.uniroot(fun, low, high))
					options(show.error.messages=FALSE)
					l[j,i]<-l[i,j]

                               } else
				{

					low <- max(low2, low3) 
					high <- min(high2, high3)					
					l[i, j] <- try(my.uniroot(fun, low, high))
					options(show.error.messages=FALSE)
					l[j,i]<-l[i,j]
				}
        
        
      
    }
    
  }
 diag(l)<- 1
 return(round(l, 3))
 }
