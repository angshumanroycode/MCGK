mcgk.test.preset<-function(n,d,parameters=NULL,B.parameters=1000){
  if(is.null(parameters))
    parameters=PARAMETERS(n-1,d,B.parameters)
  return(CGKPRESET(n-1,d,parameters))
}

mcgk.dist.test<-function(Xlist=NULL,Dlist=NULL,alpha=0.05,parameters=NULL,B.parameters=1000,
                         B=1000,test.preset=NULL){
  if(!is.null(Xlist))
    Dlist=lapply(Xlist,dist)
  Dlist=lapply(Dlist,as.matrix)
  ns=as.numeric(sapply(Dlist,dim))
  n=ns[1]
  if(sum(n!=ns)) stop("Invalid Input")
  d=length(Dlist)
  if(is.null(test.preset))
    test.preset=mcgk.test.preset(n,d,parameters,B.parameters)
  result0=MCGKTEST(Dlist,n,d,test.preset,B)
  val=result0[[1]]
  permval=result0[[2]]
  cutoff=apply(permval,1,function(x) as.numeric(quantile(x,1-alpha)))
  pval=(rowSums(permval>val)+1)/(B+1)
  result=list(val[1],val[2],cutoff[1],cutoff[2],pval[1],pval[2])
  names(result)=c("Tdist.sum.stat","Tdist.max.stat","Tdist.sum.cutoff",
                  "Tdist.max.cutoff","Tdist.sum.pvalue","Tdist.max.pvalue")
  return(result)
}

mcgk.proj.test<-function(Xlist=NULL,Ilist=NULL,alpha=0.05,parameters=NULL,B.parameters=1000,
                         B=1000,test.preset=NULL){
  
  if(!is.null(Xlist)){
    Xlist=lapply(Xlist,as.matrix)
    Ilist=lapply(Xlist,function(x) x%*%t(x))
  }
  ns=as.numeric(sapply(Ilist,dim))
  n=ns[1]
  if(sum(n!=ns)) stop("Invalid Input")
  d=length(Ilist)
  if(is.null(test.preset))
    test.preset=mcgk.test.preset(n,d,parameters,B.parameters)
  result0=MCGKTEST(Ilist,n,d,test.preset,B)
  val=result0[[1]]
  permval=result0[[2]]
  cutoff=apply(permval,1,function(x) as.numeric(quantile(x,1-alpha)))
  pval=(rowSums(permval>val)+1)/(B+1)
  result=list(val[1],val[2],cutoff[1],cutoff[2],pval[1],pval[2])
  names(result)=c("Tproj.sum.stat","Tproj.max.stat","Tproj.sum.cutoff",
                  "Tproj.max.cutoff","Tproj.sum.pvalue","Tproj.max.pvalue")
  return(result)
}