#----------------------------------------------------------------------------
require(MultiLCIRT)
#----------------------------------------------------------------------------
est_biv_LT <- function(Y,B,XX,WW,Z=NULL,k,tol=10^-8,maxit=1000,st.err=FALSE,
                                           start=0,Be=NULL,Ga=NULL,De=NULL,si2=NULL,rho=NULL,tolNR=0.01){

# estimate finite mixture bivariate lantent trajectory model
#
# INPUT:
# Y          = matrix of continuous data
# B          = matrix of binary indicator 
# XX         = matrix of covariates affecting Y
# WW         = matrix of covariates affecting B 
# Z          = matrix of covariates without intercept for the weights (optional)
# k          = number of mixture components
# r          = order ot the polynomial
# tol        = tolerance (optional)
# maxit      = max number of iterations (optional)
# st.err     = to compute standard errors (optional)
# start      = 0 for random deterministic starting, 1 for random, 2 for external (optional)

# preliminaries ----
  ncovX = dim(XX)[3]
  ncovW = dim(WW)[3]
  Y1 = Y; Y1[B==0] = NA
  sp = rowMeans(Y1,na.rm=T)
  sp = sp[!is.na(sp)]
  if(start==0){
    Be = quantile(sp,probs=(1:k)/(k+1))
    Be = rbind(Be,matrix(0,ncovX-1,k))
    mb = mean(B,na.rm=T)
    Ga = matrix(c(log(mb/(1-mb)),rep(0,ncovW-1)),ncovW,k)
    si2 = var(Y[B==1 & !is.na(B)]) 
    rho=0
  }else if(start==1){
    Be = rnorm(k)
    Be = rbind(Be,matrix(0,ncovX-1,k))
    Ga = matrix(c(-abs(rnorm((1))),rep(0,ncovW-1)),ncovW,k)
    rho = 2*runif(1)-1
    si2=runif(1)*5
  }
  n = nrow(Y)
  TT = ncol(Y)
  DD = array(0,c(2,ncovW,TT,n))
  for(i in 1:n) for(l in 1:TT) DD[2,,l,i] = WW[i,l,]
  if(k==1){
    ncov=0
  }else{
    if(is.null(Z)) Z = matrix(1,n,1) else Z = cbind(intercept=1,Z)
    ncov = ncol(Z)
  }
  if(ncov>0){
    cov_names = colnames(Z)
    out = aggr_data(Z)
    Zdis = out$data_dis; Zlabel = out$label; ndis = nrow(Zdis)
    ZZ = array(0,c(k-1,ncov*(k-1),ndis))
    I2 = diag(k-1)
    for(i in 1:ndis) ZZ[,,i] = I2%x%t(Zdis[i,])
  }
  if(ncov==0){
    de = NULL; ZZ = NULL; piv = rep(1,k)/k
  }else{
    if(start==0){
      de = NULL
      Piv = matrix(1/k,n,k)
    }
    if(start==1){
      de = NULL
      Piv = matrix(runif(k),n,k)
      Piv = Piv/rowSums(Piv)
    }
    if(start==2){
      de = as.vector(De)
      Piv = prob_multi_glob(ZZ,model="m",de,Zlabel)$P
      Be = Be
      Ga = Ga
      si2 = si2
      rho = rho
    }
  }

# compute initial log-likelihood ----
  si = sqrt(si2)
  lrho = (rho+1)/2; lrho = log(lrho/(1-lrho))
  th = c(as.vector(Be),as.vector(Ga),log(si2),lrho) 
  if(ncov==0) Piv = rep(1,n)%o%piv
  if(start==1) W = t(rmultinom(n,1,rep(1/k,k))) else W = Piv
  lPc = lk_sel_comp(th,n,k,ncovW,ncovX,DD,LL,WW,XX,W,Y,B,output=TRUE)$lPc
  mlPc = apply(lPc,1,max)
  lPc1 = lPc-mlPc
  if(ncov==0) Piv = rep(1,n)%o%piv
  Pj1 = exp(lPc1)*Piv    
  pm1 = rowSums(Pj1)
  W = Pj1/pm1
  lk = sum(log(pm1))+sum(mlPc)
  print("Initial parameter values")
  print(round(th,4))  

# iterate until convergence ----
  it = 0; lko = lk
  t0 = proc.time()
  while((abs(lk-lko)/abs(lko)>tol || it==0) && it<maxit){
    it = it+1; lko = lk
    print("EM iteration")
    print(it)
# update piv
    if(ncov==0){
      Piv = matrix(n,1)
    }else{
      out = est_multi_glob(W,ZZ,"m",Zlabel,de)
      de = out$be; Piv = out$P
    }
    piv = colMeans(Piv)
    th = optim(th, lk_sel_comp, gr = sc_sel_comp,n,k,ncovW,ncovX,DD,LL,WW,XX,W,Y,B,output=FALSE, method = "BFGS", control = list(fnscale=-1,trace=1,REPORT=1))$par
# compute log-likelihood
    lPc = lk_sel_comp(th,n,k,ncovW,ncovX,DD,LL,WW,XX,W,Y,B,output=TRUE)$lPc
    mlPc = apply(lPc,1,max)
    lPc1 = lPc-mlPc
    if(ncov==0) Piv = rep(1,n)%o%piv
    Pj1 = exp(lPc1)*Piv    
    pm1 = rowSums(Pj1)
    W = Pj1/pm1
    lk = sum(log(pm1))+sum(mlPc)
# display likelihood
    tmp = (proc.time()-t0)/it
    print(c(k,it,lk,lk-lko,tmp[3]))
    print(round(th,4))
# Newton-Raphson acceleration
    if((lk-lko)/abs(lko)<tolNR){
      print("Newton-Raphson acceleration")
#separate parameters
      Be = matrix(th[1:(k*ncovX)],ncovX,k); th = th[-(1:(k*ncovX))]  
      Ga = array(th[1:(k*ncovW)],c(ncovW,k)); th = th[-(1:(k*ncovW))]
      lsi2 = th[1]; th = th[-1];  si2 = exp(lsi2)
      lrho = th; rho = (exp(lrho)/(1+exp(lrho)))*2-1
# optimization
      th = c(as.vector(Be),as.vector(Ga),de,log(si2),lrho)
      out = lk_sel(th,n,k,ncovW,ncovX,ncov,DD,LL,WW,XX,ZZ,Zlabel,W,Y,B)
      outs = sc_sel(th,n,k,ncovW,ncovX,ncov,DD,LL,WW,XX,ZZ,Zlabel,W,Y,B,output=FALSE)
      th = optim(th, lk_sel,gr=sc_sel,n,k,ncovW,ncovX,ncov,DD,LL,WW,XX,ZZ,Zlabel,W,Y,B,method = "BFGS", control = list(fnscale=-1,trace=1,REPORT=1))$par
# separate parameters
      Be = matrix(th[1:(k*ncovX)],ncovX,k); th = th[-(1:(k*ncovX))]  
      Ga = array(th[1:(k*ncovW)],c(ncovW,k)); th = th[-(1:(k*ncovW))]  
      if(ncov == 0){
        de = NULL
      }else{
        de = th[1:(ncov*(k-1))]
        th = th[-(1:(ncov*(k-1)))]
      }
      lsi2 = th[1]; th = th[-1];  si2 = exp(lsi2)
      lrho = th; rho = (exp(lrho)/(1+exp(lrho)))*2-1
      th = c(as.vector(Be),as.vector(Ga),log(si2),lrho) 
# compute log-likelihood
      if(ncov==0) Piv = matrix(1,n) else Piv = prob_multi_glob(ZZ,model="m",de,Zlabel)$P
      lPc = lk_sel_comp(th,n,k,ncovW,ncovX,DD,LL,WW,XX,W,Y,B,output=TRUE)$lPc
      mlPc = apply(lPc,1,max)
      lPc1 = lPc-mlPc
      if(ncov==0) Piv = rep(1,n)%o%piv
      Pj1 = exp(lPc1)*Piv    
      pm1 = rowSums(Pj1)
      W = Pj1/pm1
      lk = sum(log(pm1))+sum(mlPc)
# display likelihood
      tmp = (proc.time()-t0)/it
      print(c(k,it,lk,lk-lko,tmp[3]))
      print(round(th,4))
    }
  }

# final output ----
# clustering
  clus = apply(W,1,which.max)
  table(clus)
# separate parameters
  Be = matrix(th[1:(k*ncovX)],ncovX,k); th = th[-(1:(k*ncovX))]  
  Ga = array(th[1:(k*ncovW)],c(ncovW,k)); th = th[-(1:(k*ncovW))]
  lsi2 = th[1]; th = th[-1];  si2 = exp(lsi2)
  lrho = th; rho = (exp(lrho)/(1+exp(lrho)))*2-1
  if(ncov>0) De = matrix(de,ncol=k-1)
# bic
  if(ncov==0) np = (k-1) else np = (k-1)*ncov
  np = np+k*ncovX+k*ncovX
  aic = -2*lk+2*np
  bic = -2*lk+log(n)*np
  ent = -sum(W*log(pmax(W,10^-100)))
  if(st.err){
    print("compute s.e.")
    th = c(as.vector(Be),as.vector(Ga),de,log(si2),lrho)  
    out = lk_sel(th,n,k,ncovW,ncovX,ncov,DD,LL,WW,XX,ZZ,Zlabel,W,Y,B)
    outs = sc_sel(th,n,k,ncovW,ncovX,ncov,DD,LL,WW,XX,ZZ,Zlabel,W,Y,B,output=TRUE)
    scn = NULL; Jn = NULL
    for(j in 1:length(th)){
      th1 = th; th1[j] = th1[j]+10^-6
      out1 = lk_sel(th1,n,k,ncovW,ncovX,ncov,DD,LL,WW,XX,ZZ,Zlabel,W,Y,B)
      out1s = sc_sel(th1,n,k,ncovW,ncovX,ncov,DD,LL,WW,XX,ZZ,Zlabel,W,Y,B,output=TRUE)
      scn = c(scn,(out1-out)*10^6)
      Jn = cbind(Jn,(out1s$sc-outs$sc)*10^6)
    }
    print(c(lk,out,lk-out))
    print(round(cbind(sc = outs$sc,scn=scn,diff=outs$sc-scn),4))
    Jn = -(Jn+t(Jn))/2
    iJn = ginv(Jn)
    if(any(diag(iJn)<0)) print("information matrix with negative diagonal elements")
    se = sqrt(abs(diag(iJn)))
    seBe = matrix(se[1:(k*ncovX)],ncovX,k); se = se[-(1:(k*ncovX))]
    seGa = matrix(se[1:(k*ncovW)],ncovW,k)
    se = se[-(1:(k*ncovW))]
    if(ncov>0){
      sede = se[1:(ncov*(k-1))]; se = se[-(1:(ncov*(k-1)))]
      seDe = matrix(sede,ncol=k-1)
    }
    selsi2 = se[1]; se=se[-1]
    selrho = se  
  }

# output list ----
  if(ncov==0){
    out = list(Be=Be,Ga=Ga,piv=piv,si2=si2,rho=rho,clus=clus,lk=lk,np=np,aic=aic,bic=bic,ent=ent,it=it)
  }else{ 
    out = list(Be=Be,Ga=Ga,De=De,Piv=Piv,piv=piv,si2=si2,rho=rho,clus=clus,lk=lk,np=np,aic=aic,bic=bic,ent=ent,it=it)
  }
  if(st.err){
    if(ncov==0){
      out$seBe=seBe; out$seGa=seGa; out$selsi2 = selsi2; out$selrho=selrho; out$Jn=Jn
    }else{
      out$seBe=seBe; out$seGa=seGa; out$seDe=seDe; out$selsi2 = selsi2; out$selrho=selrho; out$Jn=Jn
    }
  }
  return(out)

}