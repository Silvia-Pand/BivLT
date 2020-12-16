# *** Latent trajectory bivariate model ***

rm(list=ls())
source("est_biv_LT.R") #estimation function
source("lk_sel_comp.R")
source("sc_sel_comp.R")
source("lk_sel.R")
source("sc_sel.R")

load("ExampleData.RData") # load data

####

## Prepare data

# Response variables
Y <- Ys # matrix of continuous data (dimension n x TT)
dim(Y)
head(Y) 
B <- Bs # matrix of binary data (dimension n x TT)
dim(B)
head(B)

n <- dim(Y)[1]  #sample size
TT <- dim(Y)[2] #number of time occasions

# Covariates
XX <- XX1 # matrix of time-varying covariates (plus intercept) affecting the continuous outocome of dimension (n x TT x (ncovX+1))
dim(XX1)
WW <- WW1 # matrix of time-varying  covariates (plus intercept) affecting the binary outocome of dimension (n x TT x (ncovW+1))
dim(WW1)
ZZ <- ZZ1 # matrix of time-constant covariates (without intercept) for the weights of dimension (n x ncovZ)
dim(ZZ)

ncovX <- dim(XX)[3]-1 
ncovX
ncovW <- dim(WW)[3]-1
ncovW
ncovZ <- dim(ZZ)[2]
ncovZ

# include a second-order polynomial for year (year and year^2)
XXn <- array(0,c(n,TT,ncovX+3))
XXn[,,1:(ncovX+1)] <- XX
XXn[,,ncovX+2] <- rep(1,n)%o%((1:TT)-mean(1:TT))
XXn[,,ncovX+3] <- XXn[,,ncovX+2]^2/10
dim(XXn)
head(XXn[,1,])

WWn <- array(0,c(n,TT,ncovW+3))
WWn[,,1:(ncovW+1)] <- WW
WWn[,,ncovW+2] <- rep(1,n)%o%((1:TT)-mean(1:TT))
WWn[,,ncovW+3] <- WWn[,,ncovW+2]^2/10
dim(WWn)
head(WWn[,1,])


####

## Model estimation

k <- 3    #number of latent classes

# Model estimation with deterministic initialization
out= est_biv_LT(Y, B, XXn, WWn, Z = ZZ, k = k, start = 0,st.err=TRUE)

# Model estimation with random initialization
set.seed(321)
outr = list()
for(i in 1:5) outr[[i]] = try(est_biv_LT(Y, B, XXn, WWn, Z = ZZ, k = k, start = 1,st.err=TRUE))

save.image("ResultsSimData.RData")

####

## Display output
load("ResultsSimData.RData")

LK = rep(0,6)
LK[1]  = out$lk
for(i in 1:5) LK[i+1] = outr[[i]]$lk
LK

BIC = rep(0,6)
BIC[1]  = out$bic
for(i in 1:5) BIC[i+1] = outr[[i]]$bic
BIC
# Remark: there are not particular problems of local maxima
# We retain model with deterministic initialization (out) as the best one

# Averaged weights
colMeans(out$Piv)

out$rho # correlation between Y and B
out$si2 

## Effects of covariates
# effects on B (binary outcome)
coefGa <- round(out$Ga, 4) # dimension (ncovW x k): one coefficient for every latent class
seGa <- round(out$seGa, 4)
coefGa 
seGa


# effects on Y (continuous outcome)
coefBe = round(out$Be, 4)  # dimension (ncovX x k): one coefficient for every latent class
seBe = round(out$seBe, 4)
coefBe
seBe

# effects on latent classes
coefDe = round(out$De, 3) # dimension (ncovZ x (k-1)); one coefficient for every latent class minus 1; reference: class 1
seDe = round(out$seDe, 3)
coefDe
seDe

save.image("ResultsSimData.RData")

