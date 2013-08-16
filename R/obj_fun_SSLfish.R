obj_fun_SSLfish<-function(pars,fixed.vars){
  print(pars)
  pars[1]=exp(pars[1])*1E-15
  obj=0
  Lambda.obs=rep(1,fixed.vars$tsteps.1-51)
  I.F=rep(0,fixed.vars$n.yrs)
  I.F[(fixed.vars$tsteps.1+1):(fixed.vars$tsteps.all)]=1
  for(idem in 1:2){
    for(ieff in 1:2){
      cat(paste("idem ",idem," ieff ",ieff,"\n"))
      sim.opts=list(dem=idem,eff=ieff)
      if(idem==1)par=pars[1:5]
      else par=(pars[c(1:4,6)])
      #for(isim in 1:nreps){
        sim.result=run.1.sim(par=par, sim.opts=sim.opts, n.sites=fixed.vars$n.sites, n.yrs=fixed.vars$n.yrs,burnin=fixed.vars$tsteps.1)
        Lambda.est=sim.result$Svars$Ntot[2:fixed.vars$tsteps.1]/sim.result$Svars$Ntot[1:(fixed.vars$tsteps.1-1)]
        Lambda.est=Lambda.est[51:(fixed.vars$tsteps.1-1)]
        depletion.2000=sim.result$Svars$Ntot[fixed.vars$n.yrs-10]/250000
        obj=obj+10000*sum((Lambda.est-Lambda.obs)^2)+(sum(I.F * (colMeans(sim.result$Fvars$F.fleet)-0.3)^2)+sum(I.F * (colMeans(sim.result$Fvars$F.ssl)-0.1)^2))/fixed.vars$n.sites+30*(depletion.2000-0.17)^2
        if(depletion.2000<0.05)obj=obj+(0.05-depletion.2000)*200
      #}  
    }
  }
  cat(paste("depletion 2000 ",depletion.2000," obj fun ",obj,"\n"))
  obj
}