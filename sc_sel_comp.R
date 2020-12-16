sc_sel_comp <- function(th,n,k,ncovW,ncovX,DD,LL,WW,XX,W,Y,B,output=FALSE){

fF <- function(x) g = exp(dnorm(x,log=TRUE)-pnorm(x,log=TRUE))

# *** convert parameters ***
  Be = matrix(th[1:(k*ncovX)],ncovX,k); th = th[-(1:(k*ncovX))]
  Ga = array(th[1:(k*ncovW)],c(ncovW,k)); th = th[-(1:(k*ncovW))]
  lsi2 = th[1]; th = th[-1];  si2 = exp(lsi2)
  si2 = max(si2,10^-4)
  lrho = th; rho = (exp(lrho)/(1+exp(lrho)))*2-1
  rho = max(min(rho,0.995),-0.995)

# *** compute log-likelihood ***
  s1rho2 = sqrt(1-rho^2)
  si = sqrt(si2)
  if(output) Pc = matrix(1,n,k)
  MUB = STA = array(0,c(n,TT,k))
  for(j in 1:k){
    Mub = MM = matrix(0,n,TT)
    for(h in 1:nrow(Ga)) Mub = Mub+WW[,,h]*Ga[h,j]
    for(h in 1:nrow(Be)) MM = MM+XX[,,h]*Be[h,j]
    Sta = (Y-MM)/si
    MUB[,,j] = Mub
    STA[,,j] = Sta
    if(output){
      PP1 = pnorm(Mub)
      PP0 = 1-PP1
      Tmp = dnorm(Sta)/(si*PP1)*pnorm((rho*Sta+Mub)/s1rho2)*PP1
      for(i in 1:n){
        ind0 = which(B[i,1:TT]==0)
        ind1 = which(B[i,1:TT]==1)
        if(length(ind0)>0) Pc[i,j] = prod(PP0[i,ind0])
        if(length(ind1)>0) Pc[i,j] = Pc[i,j]*prod(Tmp[i,ind1])
      }
    }
  }
#  if(ncov==0) Piv = rep(1,n)%o%piv
  if(output) lkc = sum(W*log(Pc))

# *** compute score ***
  sc = NULL
  # regression model for Y
  FFF = array(0,c(n,TT,k))
  for(j in 1:k){
    sc0 = 0
    Sta = STA[,,j]
    FF = fF((rho*Sta+MUB[,,j])/s1rho2)
    FFF[,,j] = FF
    Tmp3 = (Sta-FF*rho/s1rho2)/si
    for(i in 1:n){
      Xii = XX[i,,]
      ind1 = which(B[i,]==1)   # ***
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
# for si2
  sc4  = 0
  dsi = 1/2*si
  for(j in 1:k){
    Sta = STA[,,j]
    Tmp2 = FFF[,,j]
    Tmp3 = Sta^2-1-rho/s1rho2*Sta*Tmp2
    for(i in 1:n){
      indl = which(B[i,]==1)# ***
      sc4 = sc4 + 1/si*W[i,j]*sum(Tmp3[i,indl])*dsi
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
      indl = which(B[i,]==1) # ***
      sc3 = sc3+ W[i,j]*(Tmp2[i,indl]%*%Tmp3[i,indl])*drho
    }
  }
  sc = c(sc,sc3)

# *** output ----
  if(any(is.nan(sc))){print("NaN sc"); print(sc)}
  if(output) out = list(lkc = lkc, sc = sc)
  else  out = sc
  return(out)

}
