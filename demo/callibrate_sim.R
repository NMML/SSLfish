require(gtools)
require(compiler)
source("c:/users/paul.conn/git/SSLfish/R/util_funcs.R")
source("c:/users/paul.conn/git/SSLfish/R/sim_funcs.R")
source("c:/users/paul.conn/git/SSLfish/R/obj_fun_SSLfish.R")
tsteps.1=150
tsteps.2=25
tsteps.3=20
tsteps.all=tsteps.1+tsteps.2+tsteps.3
n.sites=37
fixed.vars=list(tsteps.1=tsteps.1,tsteps.2=tsteps.2,tsteps.3=tsteps.3,n.sites=n.sites,tsteps.all=tsteps.all,n.yrs=tsteps.all+1)
pars=c(0.104,3.344,0.771,1.298,-0.038,-1.286) #delta, log(q), log(alpha0),log(beta0),log(-alpha1),log(-beta1)
Out<-optim(par=pars,fn=obj_fun_SSLfish,method="SANN",fixed.vars=fixed.vars)
print(Out$par)
summary(Out)
save(Out,file="sim_anneal_Out.Rdat")