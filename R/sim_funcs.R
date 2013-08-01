############################
### Simulation functions ###
############################

run.1.sim=function(par, n.sites, n.yrs){
  # par = c(p, delta, alpha0, alpha1, beta0, beta1,...)
  SSLvital=list(p=par[1], delta=par[2], alpha0=par[3], alpha1=par[4], beta0=par[5], beta1=par[6])
  Svars=set.ssl.par(n.sites,n.yrs,SSLvital)
  Fvars=set.fish.par(...)
  Fvars=???
  for(t in 1:(n.yrs-1)){
      Svars=project.ssl.1yr(Svars,t)
      Fvars=projects.fish.1yr(Fvars,Svars,n.sites,t)
      for(k in 1:n.sites){
        Svars$f[k,t+1,]=invLogit(Svars$beta0+beta1*(Svars$B[k,t]/Fvars$B[k,t]) + rnorm(32,0,1.0E-8))
        Svars$S[k,t+1,]=invLogit(Svars$alpha0+alpha1*(Svars$B[k,t]/Fvars$B[k,t]) + rnorm(32,0,1.0E-8))
      }
  }
}

##########################
### Sea lion functions ###
##########################

project.ssl.1yr=function(Svars, yr){
  # yr = current year i.e., projects from yr to yr+1
  for(k in 1:Svars$n.sites){
    Svars$N[k,yr+1,]=apply(Svars$A[k,,],1,function(a,N){sum(rbinom(64,N,a))}, N=Svars$N[k,yr,])  
    Svars$B[k,yr+1]=crossprod(unlist(Svars$mass),Svars$N[k,yr+1,])
    Svars$I.pup[k,yr+1]=rbinom(1, size=Svars$N[k,yr+1,1]+Svars$N[k,yr+1,33], prob=0.95)
    N.2p = rbinom(1, size=sum(Svars$N[k,yr+1,3:32])+sum(Svars$N[k,yr+1,35:64]), prob=Svars$p)
    Svars$I.np[k,yr+1]=max(0,round(rnorm(1, N.2p, 0.05*N.2p)))
  }  
  return(Svars)
}

# Set sea lion survival and fecundity vectors
set.ssl.pars=function(n.sites, n.yrs, SSLvital=NULL){
  attach(HFYS_appendix_C)
  if(is.null(SSLvital)){
    #data(HFYS_appendix_C) #Switch this on when package is built; remove next line
    #HFYS_appendix_C <- read.csv("data/HFYS_appendix_C.csv")
    SSLvital=list(alpha0=logit(HFYS_appendix_C$S.HFYS),
                  alpha1=-1,
                  beta0=c(logit(2*HFYS_appendix_C$f.HFYS[-1]),-Inf),
                  beta1=-1,
                  delta=1E-8,
                  p=0.5)
  }
  Svars=list(n.sites=n.sites, n.yrs=n.yrs,Ntot=round(240000*rep(1/n.sites, n.sites)), p=SSLvital$p, delta=SSLvital$delta,
    alpha0=SSLvital$alpha0, alpha1=SSLvital$alpha1, beta0=SSLvital$beta0, beta1=SSLvital$beta1, 
               mal2femSurv=c(c(HFYS_appendix_C$S.MALE.CP/HFYS_appendix_C$S.CP)[-32],0),
    mass=list(fem=c(20,richards(1:31,A=287.829,m=-0.690,S0=1.2E-04,t=4.225)),
               male=c(22,richards(1:31,A=681.112,m=8.041,S0=101.148,t=12.365))),
    S=array(0, dim=c(n.sites, n.yrs, 32)),
    f=array(0,dim=c(n.sites, n.yrs, 32)),
    N=array(0,dim=c(n.sites, n.yrs, 2*32)),
    B=matrix(0,n.sites,n.yrs),
    I.np=matrix(0,n.sites,n.yrs),
    I.pup=matrix(0,n.sites,n.yrs)
  )
  # initialize SSL survival and fecundity arrays (dim = site x year x age) 
  for(i in 1:n.sites) Svars$S[i,1,]=invLogit(Svars$alpha0)
  for(i in 1:n.sites) Svars$f[i,1,]=invLogit(Svars$beta0)
  Svars$A=set.A.array(Svars,1)
  for(i in 1:n.sites) Svars$N[i,1,]=round(Svars$Ntot[i]*Re(eigen(Svars$A[i,,])$vectors[,1]/sum(eigen(Svars$A[i,,])$vectors[,1])))
  for(i in 1:n.sites) Svars$B[i,1]=crossprod(unlist(Svars$mass),Svars$N[i,1,])
  for(i in 1:n.sites) Svars$I.pup[i,1] = rbinom(1,size=Svars$N[i,1,1]+Svars$N[i,1,33],prob=0.95)
  for(i in 1:n.sites){
    N.2p=rbinom(1, size=sum(Svars$N[i,1,3:32])+sum(Svars$N[i,1,35:64]), prob=Svars$p)
    Svars$I.np[i,1] = max(0,round(rnorm(1, N.2p, 0.05*N.2p)))
  }
  return(Svars)
}

# Create SSL Leslie projection matrix from life history params (dim = site x age x age)
set.A.array <- function(Svars,yr){
  A <- array(0,dim=c(Svars$n.sites,2*32,2*32))
  for(i in 1:Svars$n.sites){
    A[i,1,1:32] <- 0.5*Svars$f[i,yr,]*Svars$S[i,yr,]
    A[i,2:32,1:31] <- diag(Svars$S[i,yr,-32])
    A[i,33,1:32] <- 0.5*Svars$f[i,yr,]*Svars$S[i,yr,]
    A[i,34:(2*32),33:(2*32-1)] <- diag(c(Svars$mal2femSurv*Svars$S[i,yr,])[-32])
  } 
  return(A)
}


######################
### Fish functions ###
######################

get.recruitment <- function(Fvars,SSB){
  0.8*Fvars$R0*Fvars$h*SSB/(0.2*Fvars$phi0*Fvars$R0*(1-Fvars$h)+(Fvars$h-0.2)*SSB)*exp(rnorm(1,0,Fvars$sigma_R))
}

#compute unfished spawning biomass per recruit
unfished.bpr<-function(Fvars){
  n.age=length(Fvars$Mat)
  N=rep(1,n.age)
  #for(i in 1:n.age)N[i]=exp(-Fvars$M*i)
  N = exp(-Fvars$M*c(1:n.age))
  N[n.age]=N[n.age]/(1-exp(-Fvars$M))
  sum(N*Fvars$Mat*Fvars$Wgt)
}

init.fish<-function(Fvars,isite){
  N=rep(1,Fvars$n.age)
  #for(i in 1:Fvars$n.age)N[i]=exp(-Fvars$M*i)
  N = exp(-Fvars$M*c(1:n.age))
  N[n.age]=N[Fvars$n.age]/(1-exp(-Fvars$M))
  N*Fvars$R0[isite]
}

set.fish.pars<-function(n.sites,n.yrs){
  Fvars=list(a=c(1:10),M=0.2,h=0.8,R0=runif(n.sites,100000,5000000),sigma_R=0.5,mat.a=5,mat.b=1,wgt.a=5,wgt.b=0.7,eta1=1.3,eta2=1,kappa1=4,kappa2=10,eta1.SSL=1.3,eta2.SSL=1,kappa1.SSL=2,kappa2.SSL=8)
  Fvars$n.age=length(Fvars$a)
  Fvars$Mat=logistic(A=Fvars$a,a=Fvars$mat.a,b=Fvars$mat.b)
  Fvars$Wgt=logistic(A=Fvars$a,a=Fvars$wgt.a,b=Fvars$wgt.b)
  Fvars$Sel=double.logistic(A=Fvars$a,eta1=Fvars$eta1,eta2=Fvars$eta2,kappa1=Fvars$kappa1,kappa2=Fvars$kappa2)
  Fvars$Sel.SSL=double.logistic(A=Fvars$a,eta1=Fvars$eta1.SSL,eta2=Fvars$eta2.SSL,kappa1=Fvars$kappa1.SSL,kappa2=Fvars$kappa2.SSL)
  Fvars$phi0=unfished.bpr(Fvars)
  Fvars$N=array(0,dim=c(n.sites,n.yrs,Fvars$n.age))
  #initialize N at each site (array; dim = sites x years x age)
  for(isite in 1:n.sites)Fvars$N[isite,1,]=init.fish(Fvars,isite)
  #initialize Z - include natural mortality to start with
  Fvars$Z=array(Fvars$M,dim=c(n.sites,n.yrs,Fvars$n.age))
  #initialize SSB (matrix; dim = sites x years)
  Fvars$SSB=matrix(0,n.sites,n.yrs)
  Fvars$B=Fvars$SSB  #exploitable biomass
  for(isite in 1:n.sites){
    Fvars$SSB[isite,1]=sum(Fvars$N[isite,1,]*Fvars$Mat*Fvars$Wgt)
    Fvars$B[isite,1]=sum(Fvars$N[isite,1,]*Fvars$Wgt*Fvars$Sel)
  }
  Fvars
}

set.fleet.pars<-function(n.sites,n.yrs,Fvars,effort.opt){
  FLvars=list(effort=matrix(0,n.sites,n.yrs))
  if(effort.opt==1)FLvars$Effort=t(rdirichlet(n.yrs,alpha=rep(1,n.sites)))
  else FLvars$Effort[,1]=Fvars$B[,1]/sum(Fvars$B[,1])
  FLvars
}

project.fish.1yr<-function(Fvars,Svars,FLvars,n.sites,iyr){
  #add in Stellar sea lion and  mortality sources
  Fvars$Z[,iyr-1,]=Fvars$Z[,iyr-1,]+Fvars$delta*apply(Fvars$N[,iyr-1,],1,'sum')*Svars$B[iyr-1]*Fvars$Sel.SSL+FLvars$q*(FLvars$Effort[,iyr-1] %o% Fvars$Sel)
  Fvars$N[,iyr,1]=get.recruitment(Fvars,SSB=Fvars$SSB[,iyr-1])
  for(isite in 1:n.sites)Fvars$N[isite,iyr,2:Fvars$n.age]=Fvars$N[isite,iyr-1,1:(Fvars$n.age-1)]*Fvars$Z[isite,iyr-1,1:(Fvars$n.age-1)]
  Fvars$N[,iyr,Fvars$n.age]=Fvars$N[,iyr,Fvars$n.age]+Fvars$N[,iyr-1,Fvars$n.age]*Fvars$Z[,iyr-1,Fvars$n.age]  #plus group survival
  #update B, SSB, etc.
  Fvars$SSB[,iyr]=Fvars$N[,iyr,]%*%(Fvars$Wgt*Fvars$Mat)
  Fvars$B[,iyr]=Fvars$N[,iyr,]%*%(Fvars$Wgt*Fvars$Sel)
  Fvars
} 


