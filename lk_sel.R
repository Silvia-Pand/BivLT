lk_sel<- function(th,n,k,ncovW,ncovX,ncov,DD,LL,WW,XX,ZZ,Zlabel,W,Y,B){

# *** convert parameters ***
  Be = matrix(th[1:(k*ncovX)],ncovX,k); th = th[-(1:(k*ncovX))]  
  Ga = array(th[1:(k*ncovW)],c(ncovW,k)); th = th[-(1:(k*ncovW))]  
  if (ncov == 0) de = NULL
  else if(k>1){de = th[1:(ncov*(k-1))]; th = th[-(1:(ncov*(k-1)))]}
  lsi2 = th[1]; th = th[-1];  si2 = exp(lsi2)
  lrho = th; rho = (exp(lrho)/(1+exp(lrho)))*2-1

# *** compute log-likelihood ***
  if(ncov == 0){
    piv = colSums(W)/dim(W)[1]
    Piv = rep(1,dim(W)[1])%o%piv  
  }else{
    if(k==1){
      Piv = rep(1,n)
    }else{
      Piv = prob_multi_glob(ZZ,model="m",de,Zlabel)$P
    }
  }
  si = sqrt(si2)
  Pc = matrix(1,n,k)
  for(j in 1:k){
    Mub = MM = matrix(0,n,TT)
    for(h in 1:nrow(Ga)) Mub = Mub+WW[,,h]*Ga[h,j]
    for(h in 1:nrow(Be)) MM = MM+XX[,,h]*Be[h,j]
    PP1 = pnorm(Mub)
    PP0 = 1-PP1
    Sta = (Y-MM)/si
    Tmp2 = dnorm(Sta)/(si*PP1)*pnorm((rho*Sta+Mub)/sqrt(1-rho^2))*PP1
    for(i in 1:n){
      ind0 = which(B[i,]==0)
      if(length(ind0)>0) Pc[i,j] = prod(PP0[i,ind0])
      ind1 = which(B[i,]==1)
      if(length(ind1)>0) Pc[i,j] = Pc[i,j]*prod(Tmp2[i,ind1])
    }
  }
#  if(ncov==0) Piv = rep(1,n)%o%piv
  Pj = Pc*Piv
  pm = rowSums(Pj) 
  W = Pj/pm
  lk = sum(log(pm))
#  flk = -lk
  return(lk)

}