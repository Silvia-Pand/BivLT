sc_sel <- function(th,n,k,ncovW,ncovX,ncov,DD,LL,WW,XX,ZZ,Zlabel,W,Y,B,output=FALSE){

  fF <- function(x) g = dnorm(x)/pnorm(x)

# *** convert parameters ***
  Be = matrix(th[1:(k*ncovX)],ncovX,k); th = th[-(1:(k*ncovX))]
  Ga = array(th[1:(k*ncovW)],c(ncovW,k)); th = th[-(1:(k*ncovW))]
  if (ncov == 0) de = NULL
  else if(k>1){de = th[1:(ncov*(k-1))]; th = th[-(1:(ncov*(k-1)))]}
  lsi2 = th[1]; th = th[-1];  si2 = exp(lsi2)
  lrho = th; rho = (exp(lrho)/(1+exp(lrho)))*2-1
  s1rho2 = sqrt(1-rho^2)

  # *** compute log-likelihood ***
  if (ncov == 0) {
    piv = colSums(W)/dim(W)[1]
    Piv = rep(1,dim(W)[1])%o%piv  
  } else if(k==1) Piv = rep(1,n) else Piv = prob_multi_glob(ZZ,model="m",de,Zlabel)$P
  si = sqrt(si2)
  Pc = matrix(1,n,k)
  MMM = MUB = STA = array(0,c(n,TT,k))
  for(j in 1:k){
    Mub = MM = matrix(0,n,TT)
    for(h in 1:nrow(Ga)) Mub = Mub+WW[,,h]*Ga[h,j]
    for(h in 1:nrow(Be)) MM = MM+XX[,,h]*Be[h,j]
    PP1 = pnorm(Mub)
    PP0 = 1-PP1
    Sta = (Y-MM)/si
    MMM[,,j] = MM
    MUB[,,j] = Mub
    STA[,,j] = Sta
    Tmp2 = dnorm(Sta)/(si*PP1)*pnorm((rho*Sta+Mub)/s1rho2)*PP1
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

# *** compute score ***
  sc = NULL
  # regression model for Y
  FFF = array(0,c(n,TT,k))
  for(j in 1:k){
    Sta = STA[,,j]
    sc0 = 0
    FF = fF((rho*Sta+MUB[,,j])/s1rho2)
    FFF[,,j] = FF
    Tmp3 = (Sta-FF*rho/s1rho2)/si
    for(i in 1:n){    
      Xii = XX[i,,]
      ind1 = which(B[i,]==1)   # ***
      if(length(ind1)==1) sc0 = sc0 + W[i,j]*Xii[ind1,]*Tmp3[i,ind1]
      if(length(ind1)>1) sc0 = sc0+ W[i,j]*t(Xii[ind1,])%*%Tmp3[i,ind1]
    }
    sc = c(sc,sc0)
  }

# model for the probability of stopping
  for(j in 1:k){
    sc0 = sc1 = 0
    Tmp1 = dnorm(MUB[,,j],log=TRUE)-pnorm(-MUB[,,j],log.p=TRUE)
    Tmp2 = fF(MUB[,,j])
    Tmp3 = 1/s1rho2*FFF[,,j]-Tmp2
    for(i in 1:n){ 
      Wii = WW[i,,]
      ind0 = which(B[i,]==0)
      if(length(ind0)==1) sc0 = sc0 - W[i,j]*Wii[ind0,]*exp(Tmp1[i,ind0])
      if(length(ind0)>1) sc0 = sc0 - W[i,j]*t(Wii[ind0,])%*%exp(Tmp1[i,ind0])
      ind1 = which(B[i,]==1)
      if(length(ind1)==1){
        sc0 = sc0 + W[i,j]*Wii[ind1,]*Tmp2[i,ind1]
        sc1 = sc1 + W[i,j]*Wii[ind1,]*Tmp3[i,ind1]
      }
      if(length(ind1)>1){
        sc0 = sc0 + W[i,j]*t(Wii[ind1,])%*%Tmp2[i,ind1]
        sc1 = sc1 + W[i,j]*t(Wii[ind1,])%*%Tmp3[i,ind1]
      }
    }
    sc = c(sc,sc0+sc1)
  }    
# update piv
  if (ncov > 0) if(k>1) sc = c(sc,est_multi_glob(W,ZZ,"m",Zlabel,de,only_sc=T)$sc)
# for si2
  sc4 = 0
  dsi = 1/2*si
  for(j in 1:k){
    Sta = STA[,,j]
    Tmp2 = FFF[,,j]
    Tmp3 = Sta^2-1-rho/s1rho2*Sta*Tmp2
    for(i in 1:n){
      ind1 = which(B[i,]==1) # ***
      sc4 = sc4 + 1/si*W[i,j]*sum(Tmp3[i,ind1])*dsi
      #      browser()
#      sc4 = sc4 + W[i,j]*1/si*((((Y[i,indl]-Mui[indl])/si)^2-1-rho/sqrt(1-rho^2)*(Y[i,indl]-Mui[indl])/si)%*%tmp2)*dsi
    }
  }
  sc = c(sc,sc4)  

#for rho
drho = 2*exp(lrho)/(1+exp(lrho))^2
sc3 = 0
  for(j in 1:k){
    Sta = STA[,,j]
    Tmp2 = FFF[,,j]
    Tmp3 = (Sta+rho*MUB[,,j])/((1-rho^2)*s1rho2)
    for(i in 1:n){
      ind1 = which(B[i,]==1) # ***
      sc3 = sc3+ W[i,j]*(Tmp2[i,ind1]%*%Tmp3[i,ind1])*drho
    }
  }
  sc = c(sc,sc3)  
#  dflk = -sc  
#  return(sc)
  
# *** output ***
  if(output) out = list(lk=lk,sc=sc)
  else out = sc
  return(out)
}
