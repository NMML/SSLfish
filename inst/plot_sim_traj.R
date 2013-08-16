#simdriver.R - Script to conduct SSL/fish/fisheries simulations
require(gtools)
require(compiler)
require(glmmLDTS)
source('c:/users/paul.conn/git/SSLfish/R/util_funcs.R')
source('c:/users/paul.conn/git/SSLfish/R/sim_funcs.R')
setwd('c:/users/paul.conn/git/SSLfish')
every=nreps/10
tsteps.1=150  # number of years with no fishing (burn in)
tsteps.2=25   # number of years with fishing but no survey
tsteps.3=20   # number of years with fishing and survey
n.sites=37
n.sim=10
tsteps.all=tsteps.1+tsteps.2+tsteps.3
n.yrs=tsteps.all+1
dem.opt=c("fec","surv")
effort.opt=c("rand","prop")
#pars=c(1E-15,log(10),log(1),log(1),-log(1),-log(1)) #delta, q, alpha0,beta0,alpha1,beta1
pars=c(1E-15*exp(0.104),3.344,0.771,1.298,-0.038,-1.286) #delta, log(q), log(alpha0),log(beta0),log(-alpha1),log(-beta1)

idem=1
ieff=1
dem=dem.opt[idem]
eff=effort.opt[ieff]
cur.file=paste("SSLfish_",dem,"_",eff,".csv",sep='')
sim.opts=list(dem=idem,eff=ieff)
par=pars[1:5]
N=matrix(0,n.sim,n.yrs)
F.fleet=N
F.ssl=N

for(isim in 1:n.sim){
  sim.result=run.1.sim(par=par, sim.opts=sim.opts, n.sites=n.sites, n.yrs=n.yrs,burnin=tsteps.1)
  N[isim,]=sim.result$Svars$Ntot
  F.fleet[isim,]=apply(sim.result$Fvars$F.fleet,2,'mean')
  F.ssl[isim,]=apply(sim.result$Fvars$F.ssl,2,'mean')
}

plot(N[1,],ylim=c(0,250000),type="l",lwd=1.2,cex.axis=1.3,cex.lab=1.3,xlab="Year",ylab="SSL (Numbers)")
for(isim in 2:n.sim)lines(N[isim,],lwd=1.2)
lines(apply(N,2,'mean'),lwd=4)

idem=2
ieff=1
dem=dem.opt[idem]
par=(pars[c(1:4,6)])
sim.opts=list(dem=idem,eff=ieff)
for(isim in 1:n.sim){
  sim.result=run.1.sim(par=par, sim.opts=sim.opts, n.sites=n.sites, n.yrs=n.yrs,burnin=tsteps.1)
  N[isim,]=sim.result$Svars$Ntot
  F.fleet[isim,]=apply(sim.result$Fvars$F.fleet,2,'mean')
  F.ssl[isim,]=apply(sim.result$Fvars$F.ssl,2,'mean')
}

for(isim in 1:n.sim)lines(N[isim,],lwd=1.2,col='gray')
lines(apply(N,2,'mean'),lwd=4,col='gray')

legend(25,100000,c("Survival","Fecundity"),lwd=c(3,3),lty=c(1,1),cex=1.3,col=c("black","gray"))

