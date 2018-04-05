assoc_comb <-
function(x, rZiZj,rZiZjR,rZiZjY,rZiRZjR,rZiRZjY,rZiYZjY){

R<-matrix(0,nrow=3*ncol(x),ncol=3*ncol(x))

sub1<-seq(from=1, to=3*ncol(x)-2, by=3)
sub2<-seq(from=2, to=3*ncol(x)-1, by=3)
sub3<-seq(from=3, to=3*ncol(x), by=3)

R[sub1,sub1]<-rZiZj
R[sub1,sub2]<-rZiZjR
R[sub1,sub3]<-rZiZjY
R[sub2,sub1]<-t(rZiZjR)
R[sub2,sub2]<-rZiRZjR
R[sub2,sub3]<-rZiRZjY
R[sub3,sub1]<-t(rZiZjY)
R[sub3,sub2]<-t(rZiRZjY)
R[sub3,sub3]<-rZiYZjY
return(R)
}
