

calCri1=function(A,X,K)
 {
     n=dim(A)[1]
     M=matrix(0,K,K)
     n=dim(A)[1]
     T=0
   for (k in 1:K)
    for (l in 1:K)
     { ind1=which(X==k)
       ind2=which(X==l)
       M[k,l]=sum(A[ind1,ind2])
     }
   for (k in 1:K)
    for (l in 1:K)
     {
       if (M[k,l]>0)
        {
          T=T+log(M[k,l]/(sum(M[k,])*sum(M[,l])))*M[k,l]
        }
     }
   return(list(value=T,M1=M))
 }
calCri2=function(A,X,K,I,new,M)
 { 
   old=X[I]
   for (k in 1:K)
    if ((k!=new)&(k!=old))
     {
      s=sum(A[I,(X==k)]) 
      M[k,new]=M[k,new]+s
      M[new,k]=M[k,new]
      M[k,old]=M[k,old]-s
      M[old,k]=M[k,old]
     }
   s=sum(A[I,(X==old)])-A[I,I]
   t=sum(A[I,(X==new)])
   M[old,new]=M[old,new]-t+s
   M[new,old]=M[old,new]
   M[old,old]=M[old,old]-2*s-A[I,I]
   M[new,new]=M[new,new]+2*t+A[I,I]
   X[I]=new
   T=0

   for (k in 1:K)
    for (l in 1:K)
     {
       if (M[k,l]>0)
        {
          T=T+log(M[k,l]/(sum(M[k,])*sum(M[,l])))*M[k,l]
        }
     }
   return(list(value=T,M1=M))
 }
mutiExp=function(A,init,K) ##K number of clusters
 {

   n=dim(A)[1]
   LTenure=min(20,ceiling(n/4))
  
   m=5000*n ##times of iteration
   X=init  ##initial VOptimal
   XOptimal=X
   L=matrix(0,n,K)
   result=calCri1(A,X,K)
   VOptimal=result$value
   M=result$M1
   localTime=0  ##time reaching a better solution
   t=0
   while ((t<m)&(localTime<1000)) 
    {
      VLocal=-Inf
      flagMove=0
      for (i in 1:n)
       { VK=-Inf
         change=0
       for (k in 1:K)
        if ((L[i,k]==0)&(X[i]!=k)&(sum(X==X[i])>2))
         {
            flagMove=1
            result=calCri2(A,X,K,i,k,M)
            V=result$value
            if (V>VK) { VK=V 
                       change=k }
         }
       if (VK>VOptimal) {I=i
                         iterChange=change
                         VOptimal=VK
                         XOptimal=X
                         XOptimal[I]=iterChange
                         localTime=0
                         break }
       if (VK>VLocal) {I=i
                       iterChange=change
                       VLocal=VK}
       }
      if (flagMove==0) break
      localTime=localTime+1  
      result=calCri2(A,X,K,I,iterChange,M) 
      X[I]=iterChange
      M=result$M1
      L=L-1
      L[L<0]=0
      L[I,iterChange]=LTenure
    }
  return(list(XOptimal=XOptimal,VOptimal=VOptimal))
 }

####################
##Tabu Search from different starting points and different point order
####################
tabuDiffStart=function(A,K)
 {
   n=dim(A)[1]
   optimalSol=-Inf
   optimalCut=rep(0,n)
   cut=rep(0,n)

   for (j in 1:5)
    {
       rndInd=sample(1:n)
       A1=A[rndInd,rndInd]
       m=n%/%K
       init=rep(1,n)
       init1=rep(1:K, each =m)
       init[1:(m*K)]=init1
       init=sample(init)
       result=mutiExp(A1,init,K)
       cut1=result$XOptimal
       if (result$VOptimal>optimalSol)
       {  for (i in 1:n)
            { optimalCut[rndInd[i]]=cut1[i]}
          optimalSol=result$VOptimal
       }
  
    }
   return(list(optimalCut=optimalCut,optimalSol=optimalSol))
 } 
