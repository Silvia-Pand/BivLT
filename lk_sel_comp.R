lk_sel_comp<- function(th,n,k,ncovW,ncovX,DD,LL,WW,XX,W,Y,B,output=FALSE){
  
# *** convert parameters ***
  # print("theta")
  # print(th)
  Be = matrix(th[1:(k*ncovX)],ncovX,k); th = th[-(1:(k*ncovX))]
  Ga = array(th[1:(k*ncovW)],c(ncovW,k)); th = th[-(1:(k*ncovW))]
  lsi2 = th[1]; th = th[-1]; si2 = exp(lsi2)
  si2 = max(si2,10^-4)  
  lrho = th; rho = (exp(lrho)/(1+exp(lrho)))*2-1
  rho = max(min(rho,0.995),-0.995)
# *** compute log-likelihood ***
  si = sqrt(si2)
  ind0 = which(B==0)
  ind1 = which(B==1)
  lPc = matrix(0,n,k)
  for(j in 1:k){
    Mub = MM = Tmp3 = matrix(0,n,TT)
    for(h in 1:nrow(Ga)) Mub = Mub+WW[,,h]*Ga[h,j]
    for(h in 1:nrow(Be)) MM = MM+XX[,,h]*Be[h,j]
    lPP1 = pnorm(Mub,log=TRUE)
    lPP0 = pnorm(-Mub,log=TRUE)
    Tmp = (Y-MM)/si
    Tmp2 = dnorm(Tmp,log=TRUE)-(log(si)+lPP1)+pnorm((rho*Tmp+Mub)/sqrt(1-rho^2),log=TRUE)+lPP1
    Tmp3[ind0] = lPP0[ind0]
    Tmp3[ind1] = Tmp2[ind1]
    lPc[,j] = rowSums(Tmp3,na.rm=TRUE)
  }
#  if(ncov==0) Piv = rep(1,n)%o%piv
  lkc = sum(W*lPc)
  if(is.nan(lkc) | is.infinite(lkc)){ print("NaN or Inf lkc"); print(lkc)}

  # print("lkc")
  # print(lkc)
#  flk = -lk
  if(output) out = list(lkc = lkc, lPc = lPc)
  else  out = lkc
  return(out)    
}
