############################
### Simulation functions ###
############################

#' Execute one simulation run
#' 
#' This function takes the arguments of \code{par} and runs the population simulation 1 time for the purposes
#' of calibration. 
#' 
#' @param par Vector of simulation variables subject to calibration
#' @param sim.opts ???
#' @param n.sites Number of sites within the simulation
#' @param n.yrs Number of total years in the simulation
#' @param burnin ??? 
#' @details This function runs a single simulation run so that the free variables in the par vector can be calibrated
#' to provide realistic population growth and interaction between sea lions, prey and fisheries components.
#' @author Paul Conn
#' @export
#' 
run.1.sim=function(par, sim.opts, n.sites, n.yrs, burnin){
  #set.seed(12345)
  #n.sites=10
  #n.yrs=100
  #burnin=25
  #par = c(1E-15,.000001,1,1,-1)
  #sim.opts=list(dem=1,eff=1)
  if(sim.opts$dem==1)SSLvital=list(p=0.5, delta=par[1],alpha0=exp(par[3]),beta0=exp(par[4]),alpha1=-exp(par[5]), beta1=0) #put log link on alpha0, beta0 to prevent '0' Holmes estimates from becoming infinite on the log- or logit- scale
  else SSLvital=list(p=0.5, delta=par[1],alpha0=exp(par[3]),beta0=exp(par[4]),alpha1=0,beta1=-exp(par[5]))
  Svars=set.ssl.pars(n.sites,n.yrs,SSLvital)
  Fvars=set.fish.pars(n.sites,n.yrs)
  FLvars=set.fleet.pars(n.sites,n.yrs)
  if(burnin==0)FLvars=update.fleet.pars(n.sites,1,Fvars,sim.opts,burnin)
  FLvars$q=exp(par[2])
  for(t in 1:(n.yrs-1)){
      Svars=project.ssl.1yr(Svars,t) #this function projects from t to t+1
      Fvars=project.fish.1yr(Fvars,Svars,FLvars,n.sites,t+1) #this function projects from t-1 to t
      FLvars=update.fleet.pars(n.sites,t+1,Fvars,FLvars,sim.opts,burnin)
      for(k in 1:n.sites){
        Svars$f[k,t+1,]=invLogit(Svars$beta0+Svars$beta1*(Svars$B[k,t]/Fvars$X[k,t]) + rnorm(32,0,1.0E-8))
        Svars$S[k,t+1,]=invLogit(Svars$alpha0+Svars$alpha1*(Svars$B[k,t]/Fvars$X[k,t]) + rnorm(32,0,1.0E-8))
      }
      Svars$A=set.A.array(Svars,t+1)
  }
  Svars$Ntot=apply(Svars$N,2,'sum')
  Out=list(Fvars=Fvars,Svars=Svars,FLvars=FLvars)
}

##########################
### Sea lion functions ###
##########################

#' Project sea lion populations
#' 
#' This function updates the sea lion population by one year.
#' 
#' @param Svars Named list of sea lion variables craeted from \code{\link{set.ssl.pars}}.
#' @param yr Current year.
#' 
#' @details This function takes survival and natality information from the \code{Svars} object to update the 
#' sea lion population from year 'yr' to 'yr+1'. The simulated survey sample of the population is accomplished 
#' at the same time.
#' @author Devin Johnson
#' @export

project.ssl.1yr=function(Svars, yr){
  # yr = current year i.e., projects from yr to yr+1
  for(k in 1:Svars$n.sites){
    Svars$N[k,yr+1,]=apply(Svars$A[k,,],1,function(a,N){sum(rbinom(64,N,a))}, N=Svars$N[k,yr,])  
    #cur.lam=Svars$N[k,yr,]%*%Svars$A[k,1,]
    #Svars$N[k,yr+1,c(1,33)]=rpois(2,cur.lam)
    Svars$B[k,yr+1]=crossprod(unlist(Svars$mass),Svars$N[k,yr+1,])
    #Svars$I.pup[k,yr+1]=rbinom(1, size=Svars$N[k,yr+1,1]+Svars$N[k,yr+1,33], prob=0.95)
    N.p=Svars$N[k,yr+1,1]+Svars$N[k,yr+1,33]
    Svars$I.pup[k,yr+1]=rbinom(1,N.p,0.95)
    #N.2p = rbinom(1, size=sum(Svars$N[k,yr+1,3:32])+sum(Svars$N[k,yr+1,35:64]), prob=Svars$p)
    N.2p = Svars$p*(sum(Svars$N[k,yr+1,2:32])+sum(Svars$N[k,yr+1,34:64]))
    Svars$I.np[k,yr+1]=N.2p*exp(rnorm(1,0,0.05))
  }  
  return(Svars)
}


#' Set sea lion simulation variables
#' 
#' @param n.sites Number of sea lion populations
#' @param n.yrs Number of years to run simulation.
#' @param SSLvital Optional list of scaling parameters to adjust baseline survival and fecundity
#' from the HFYS 1970s values.
#' @details This function sets the initial simulation variables for the sea lion populations. As a baseline 
#' the 1970s survival and natality parameters of Holmes et al. (2007) are used to initialize the female populations.
#' These values can be adjusted via the \code{SSLvital} argument which allows the user to provide multipliers for 
#' these parameters on the logit scale. The male survival parameters are taken from 
#' Calkins and Pitcher (???). The age structured mass is calculated from the fitted Richards growth models
#' presented in Calkins and Pitcher (???) as well.
#' @author Devin Johnson and Paul Conn
#' @export
set.ssl.pars=function(n.sites, n.yrs, SSLvital=NULL){
  data(HFYS_appendix_C) #Switch this on when package is built; remove next line
  #HFYS_appendix_C <- read.csv("data/HFYS_appendix_C.csv")
  if(is.null(SSLvital)){  ##right now, this would set the scaling parameters alpha0,beta0 to 1
    SSLvital=list(alpha0=logit(HFYS_appendix_C$S.HFYS),
                  alpha1=-1,
                  beta0=c(logit(HFYS_appendix_C$f.HFYS[-1]),-Inf),
                  beta1=-1,
                  delta=1E-14,
                  p=0.5)
  }
  else{
    SSLvital$alpha0=logit(HFYS_appendix_C$S.HFYS)*SSLvital$alpha0
    SSLvital$beta0=c(logit(HFYS_appendix_C$f.HFYS[-1])*SSLvital$beta0,-Inf)             
  }
  Svars=list(n.sites=n.sites, n.yrs=n.yrs,Ntot=round(250000*rep(1/n.sites, n.sites)), p=SSLvital$p, delta=SSLvital$delta,
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
    #N.2p=rbinom(1, size=sum(Svars$N[i,1,3:32])+sum(Svars$N[i,1,35:64]), prob=Svars$p)
    N.2p=(sum(Svars$N[i,1,2:32])+sum(Svars$N[i,1,34:64]))*Svars$p
    Svars$I.np[i,1] = N.2p*exp(rnorm(1,0,0.05))
  }
  return(Svars)
}

#' Set sea lion projection matrix
#' 
#' Takes survival and natality information from the sea lion simulation variables object and create 
#' a Leslie projection matrix for each site.
#' 
#' @param Svars Named list created with the \code{\link{set.ssl.pars}} function.
#' @param yr Current year.
#' 
#' @details Takes the survival and natality information for each site in year \code{yr} and places the information into a 
#' Leslie matrix for population projection to year \code{yr}+1.
#' @author Devin Johnson
#' @export
#' 

set.A.array <- function(Svars,yr){
  A <- array(0,dim=c(Svars$n.sites,2*32,2*32))
  for(i in 1:Svars$n.sites){
    A[i,1,1:32] <- Svars$f[i,yr,]*0.949*Svars$S[i,yr,] #0.949 is neonate survival
     A[i,2:32,1:31] <- diag(Svars$S[i,yr,-32])
    A[i,33,1:32] <- Svars$f[i,yr,]*0.949*Svars$S[i,yr,]
    A[i,34:(2*32),33:(2*32-1)] <- diag(c(Svars$mal2femSurv*Svars$S[i,yr,])[-32])
  } 
  return(A)
}


######################
### Fish functions ###
######################

#' @title Calculate fish recruitment
#' @param Fvars Named list of fish simulation variables created by call to \code{\link{set.fish.pars}}
#' @param SSB Spawning biomass
#' @author Paul Conn
#' @export
#' 
get.recruitment <- function(Fvars,SSB){
  0.8*Fvars$R0*Fvars$h*SSB/(0.2*Fvars$phi0*Fvars$R0*(1-Fvars$h)+(Fvars$h-0.2)*SSB)*exp(rnorm(1,0,Fvars$sigma_R))
}

#' @title Compute unfished spawning biomass per recruit
#' @param Fvars Named list of fish simulation variables created by call to \code{\link{set.fish.pars}}
#' @author Paul Conn
#' @export
unfished.bpr<-function(Fvars){
  n.age=length(Fvars$Mat)
  N=rep(1,Fvars$n.age)
  #for(i in 1:n.age)N[i]=exp(-Fvars$M*i)
  N = exp(-Fvars$M*c(1:Fvars$n.age))
  N[Fvars$n.age]=N[Fvars$n.age]/(1-exp(-Fvars$M))
  sum(N*Fvars$Mat*Fvars$Wgt)
}

#' @title Initialize fish abundance
#' @param Fvars Named list of fish simulation variables created by call to \code{\link{set.fish.pars}}
#' @param isite Integer site number
#' @author Paul Conn
#' @export
#' 
init.fish<-function(Fvars,isite){
  N=rep(1,Fvars$n.age)
  #for(i in 1:Fvars$n.age)N[i]=exp(-Fvars$M*i)
  N = exp(-Fvars$M*c(1:Fvars$n.age))
  N[Fvars$n.age]=N[Fvars$n.age]/(1-exp(-Fvars$M))
  N*Fvars$R0[isite]
}

#' @title Initialize fish simulation variables
#' @param n.sites Number of sites in the simulation
#' @param n.yrs Number of total years in the simulation 
#' @details Add more details here...
#' @author Paul Conn
#' @export
#' 
set.fish.pars<-function(n.sites,n.yrs){
  Fvars=list(a=c(1:10),M=0.2,h=0.8,R0=runif(n.sites,1000000,5000000),sigma_R=0.5,mat.a=5,mat.b=1,wgt.a=5,wgt.b=0.7,eta1=1.3,eta2=1,kappa1=4,kappa2=10,eta1.SSL=1.3,eta2.SSL=1,kappa1.SSL=2,kappa2.SSL=8)
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
  Fvars$F.ssl=matrix(0,n.sites,n.yrs)
  Fvars$F.fleet=matrix(0,n.sites,n.yrs)
  #initialize SSB (matrix; dim = sites x years)
  Fvars$SSB=matrix(0,n.sites,n.yrs)
  Fvars$B=Fvars$SSB  #exploitable biomass
  Fvars$X=Fvars$SSB  #exploitable biomass
  for(isite in 1:n.sites){
    Fvars$SSB[isite,1]=sum(Fvars$N[isite,1,]*Fvars$Mat*Fvars$Wgt)
    Fvars$B[isite,1]=sum(Fvars$N[isite,1,]*Fvars$Wgt*Fvars$Sel)
    Fvars$X[isite,1]=sum(Fvars$N[isite,1,]*Fvars$Wgt*Fvars$Sel.SSL)
  }
  Fvars$Catch.num=matrix(0,n.sites,n.yrs)
  Fvars$Catch.wgt=Fvars$Catch.num
  Fvars
}

#' @title Initialize fishing fleet simulation variables
#' @param n.sites Number of sites in the simulation
#' @param n.yrs Number of total years in the simulation 
#' @param burnin Number of years in the simulation before 'fishing' begins.
#' @details Add more details here...
#' @author Paul Conn
#' @export
#' 
set.fleet.pars<-function(n.sites,n.yrs,burnin){
  FLvars=list(Effort=matrix(0,n.sites,n.yrs))
  FLvars
}

#' @title Update fishing fleet simulation variables
#' @param n.sites Number of sites in the simulation
#' @param iyr Current year in the simulation
#' @param Fvars Fish population variables list created by call to \code{\link{set.fish.pars}}
#' @param FLvars Fishing fleets variables list created by call to \code{\link{set.fleet.pars}}
#' @param sim.opts ???
#' @param burnin Number of years in the simulation before 'fishing' begins.
#' @author Paul Conn
#' @export
#' @import gtools
#' 
update.fleet.pars<-function(n.sites,iyr,Fvars,FLvars,sim.opts,burnin){
  if(iyr>burnin){
    if(sim.opts$eff==1)FLvars$Effort[,iyr]=rdirichlet(1,alpha=rep(1,n.sites))
    else FLvars$Effort[,iyr]=Fvars$B[,iyr]/sum(Fvars$B[,iyr])
  }
  FLvars
}

#' @title Project fish population forward 1 year
#' @param Fvars Fish population variables list created by call to \code{\link{set.fish.pars}}
#' @param Svars Sea lion population variables list created by call to \code{\link{set.ssl.pars}}
#' @param FLvars Fishing fleets variables list created by call to \code{\link{set.fleet.pars}}
#' @param n.sites Number of sites in the simulation study
#' @param iyr Current year in the simulation
#' @author Paul Conn
#' @export
#' 
project.fish.1yr<-function(Fvars,Svars,FLvars,n.sites,iyr){
  #add in Stellar sea lion and  mortality sources
  Fvars$F.fleet[,iyr-1]=FLvars$q*FLvars$Effort[,iyr-1]
  Cur.F=Fvars$F.fleet[,iyr-1] %o% Fvars$Sel
  Fvars$F.ssl[,iyr-1]=Svars$delta*(apply(Fvars$N[,iyr-1,],1,'sum')*Svars$B[,iyr-1])
  Fvars$Z[,iyr-1,]=Fvars$Z[,iyr-1,]+Fvars$F.ssl[,iyr-1]%o%Fvars$Sel.SSL+Cur.F
  Fvars$N[,iyr,1]=get.recruitment(Fvars,SSB=Fvars$SSB[,iyr-1])
  for(isite in 1:n.sites)Fvars$N[isite,iyr,2:Fvars$n.age]=Fvars$N[isite,iyr-1,1:(Fvars$n.age-1)]*exp(-Fvars$Z[isite,iyr-1,1:(Fvars$n.age-1)])
  Fvars$N[,iyr,Fvars$n.age]=Fvars$N[,iyr,Fvars$n.age]+Fvars$N[,iyr-1,Fvars$n.age]*exp(-Fvars$Z[,iyr-1,Fvars$n.age])  #plus group survival
  #update B, SSB, etc.
  Fvars$SSB[,iyr]=Fvars$N[,iyr,]%*%(Fvars$Wgt*Fvars$Mat)
  Fvars$B[,iyr]=Fvars$N[,iyr,]%*%(Fvars$Wgt*Fvars$Sel)
  Fvars$X[,iyr]=Fvars$N[,iyr,]%*%(Fvars$Wgt*Fvars$Sel.SSL)  
  Catch.age=Cur.F/Fvars$Z[,iyr-1,]*(1-exp(-Fvars$Z[,iyr-1,]))*Fvars$N[,iyr-1,]
  Fvars$Catch.num[,iyr-1]=apply(Catch.age,1,'sum')
  Fvars$Catch.wgt[,iyr-1]=Catch.age%*%Fvars$Wgt
  Fvars
} 





