
###
### Functions to increment populations by one year
###

### Sea lion functions
library(compiler)

inc.N.ssl <- function(A.list, N.ssl.list){
  
}


get.recruitment <- function(Fvars,SSB){
  0.8*Fvars$R0*Fvars$h*SSB/(0.2*Fvars$phi0*Fvars$R0*(1-Fvars$h)+(Fvars$h-0.2)*SSB)*exp(rnorm(1,0,Fvars$sigma_R))
}

#compute unfished spawning biomass per recruit
unfished.bpr<-function(Fvars){
  n.age=length(Fvars$Mat)
  N=rep(1,n.age)
  for(i in 1:n.age)N[i]=exp(-Fvars$M*i)
  N[n.age]=N[n.age]/(1-exp(-Fvars$M))
  sum(N*Fvars$Mat*Fvars$Wgt)
}

init.fish<-function(Fvars,isite){
  N=rep(1,Fvars$n.age)
  for(i in 1:Fvars$n.age)N[i]=exp(-Fvars$M*i)
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
=======
inc.S <- function(B.star.list, B.list, param.S.list){
  
}

inc.f <- function(B.star.list, B.list, param.f.list){
  
}

make.A.list <- cmpfun(function(f.list, S.list, mal2femRatio){
  K <- length(f.list)
  n.age <- length(f.list[[1]])
  out <- replicate(K, matrix(0,2*n.age,2*n.age), simplify=FALSE)
  for(i in 1:K){
    out[[i]][1,1:n.age] <- 0.5*f.list[[i]]*S.list[[i]]
    out[[i]][2:n.age, 1:(n.age-1)] <- diag(S.list[[i]][-n.age])
    out[[i]][(n.age+1),1:n.age] <- 0.5*f.list[[i]]*S.list[[i]]
    out[[i]][(n.age+2):(2*n.age), (n.age+1):(2*n.age-1)] <- diag((mal2femRatio*S.list[[i]])[-n.age])
  } 
  return(out)
})

