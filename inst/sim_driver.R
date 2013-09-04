#simdriver.R - Script to conduct SSL/fish/fisheries simulations
require(gtools)
require(compiler)
require(glmmLDTS)
source('c:/users/paul.conn/git/SSLfish/R/util_funcs.R')
source('c:/users/paul.conn/git/SSLfish/R/sim_funcs.R')
setwd('c:/users/paul.conn/git/SSLfish')
set.seed(12345)
nreps=1000 #make sure divisible by 10
every=nreps/10
tsteps.1=150  # number of years with no fishing (burn in)
tsteps.2=25   # number of years with fishing but no survey
tsteps.3=20   # number of years with fishing and survey
n.sites=37
tsteps.all=tsteps.1+tsteps.2+tsteps.3
n.yrs=tsteps.all+1
dem.opt=c("surv","fec")
effort.opt=c("rand","prop")
Result=data.frame(matrix(0,nreps,18))
colnames(Result)=c("dem","effort","p.adcount.cpue","p.adcount.catch","p.adcount.effort","p.adcount.prey","p.pupcount.cpue","p.pupcount.catch","p.pupcount.effort","p.pupcount.prey","p.adlam.cpue","p.adlam.catch","p.adlam.effort","p.adlam.prey","p.puplam.cpue","p.puplam.catch","p.puplam.effort","p.puplam.prey")
#pars=c(1E-15,log(10),log(1),log(1),-log(1),-log(1)) #delta, q, alpha0,beta0,alpha1,beta1
#pars=c(1E-15*exp(0.10),3.34,0.77,1.29,-0.04,-1.29) #delta, log(q), log(alpha0),log(beta0),log(-alpha1),log(-beta1)
pars=matrix(0,4,5)
load('sim_anneal_Out_12.Rdat')
pars[1,]=Out$par
pars[2,]=Out$par
load('sim_anneal_Out_21.Rdat')
pars[3,]=Out$par
pars[4,]=Out$par
pars[,1]=exp(pars[,1])*1E-15
#pars[1,]=c(1.7432553,2.1368623,2.0029217,0.5221643,3.1616766)
#obj,~,138000
#pars[2,]=c(2.1974836,2.5055605,4.3300752,0.3688895,4.5284402)
#obj,~158000
#pars[3,]=c(1.7692918,2.5519607,0.9883099,-0.2321302,1.2464747)
#obj,~237000
#pars[4,]=c(-0.3587023,2.3022877,1.5209649,0.3378606,1.5414819)
for(idem in 1:1){
  for(ieff in 1:1){
    dem=dem.opt[idem]
    eff=effort.opt[ieff]
    cur.file=paste("SSLfish_",dem,"_",eff,".csv",sep='')
    sim.opts=list(dem=idem,eff=ieff)
    #if(idem==1)par=pars[1:5]
    #else par=(pars[c(1:4,6)])
    par=pars[(idem-1)*2+ieff,]
    for(isim in 1:nreps){
      cat(paste("sim ",isim,"\n"))
      sim.result=run.1.sim(par=par, sim.opts=sim.opts, n.sites=n.sites, n.yrs=n.yrs,burnin=tsteps.1)
      CPUE=matrix(exp(rnorm(tsteps.3*n.sites,log(sim.result$Fvars$Catch.wgt[,(tsteps.1+tsteps.2+1):tsteps.all]/sim.result$FLvars$Effort[,(tsteps.1+tsteps.2+1):tsteps.all]),sqrt(log(0.2^2+1)))),n.sites,tsteps.3)      
      Catch=sim.result$Fvars$Catch.wgt[,(tsteps.1+tsteps.2+1):tsteps.all]
      Effort=sim.result$FLvars$Effort[,(tsteps.1+tsteps.2+1):tsteps.all]
      Prey=sim.result$Fvars$B[,(tsteps.1+tsteps.2+1):tsteps.all]
      Pups=sim.result$Svars$I.pup[,(tsteps.1+tsteps.2+1):n.yrs]
      Adults=sim.result$Svars$I.np[,(tsteps.1+tsteps.2+1):n.yrs]
      Dep.vars=c("Adults","Pups")
      Ind.vars=c("CPUE","Catch","Effort","Prey")
      Result[isim,1]=dem
      Result[isim,2]=eff
      #count analyses
      Count.df=data.frame(cbind(rep(1:n.sites,each=tsteps.3),rep(1:tsteps.3,n.sites),as.vector(t(Adults[,-1])),as.vector(t(Pups[,-1])),as.vector(t(CPUE))),as.vector(t(Catch)),as.vector(t(Effort)),as.vector(t(Prey)))
      colnames(Count.df)=c("Site","Time","Adults","Pups","CPUE","Catch","Effort","Prey")
      Count.df[,"CPUE"]=Count.df[,"CPUE"]/mean(Count.df[,"CPUE"])
      Count.df[,"Catch"]=Count.df[,"Catch"]/mean(Count.df[,"Catch"])
      Count.df[,"Effort"]=Count.df[,"Effort"]/mean(Count.df[,"Effort"])
      Count.df[,"Prey"]=Count.df[,"Prey"]/mean(Count.df[,"Prey"])
      for(idep in 1:2){
        for(iind in 1:4){
          dep=Dep.vars[idep]
          ind=Ind.vars[iind]
          fixed.formula=formula(eval(parse(text=paste(dep,'~',ind))))
          random.formula=formula(eval(parse(text=paste(dep,'~ Site'))))          
          cur.mod<-glmmLDTS(fixed.formula=fixed.formula,random.formula=random.formula,data=Count.df,timecol="Time",group.vec="Site",distribution="Poisson")          
          Result[isim,2+(idep-1)*4+iind]=sign(cur.mod$fixed.effects[2,"estimate"])*cur.mod$fixed.effects[2,"prob.t"]
          if(Result[isim,2+(idep-1)*4+iind]==0 & sign(cur.mod$fixed.effects[2,"estimate"])==-1)Result[isim,2+(idep-1)*4+iind]=-2
        }
      }
      #lambda analyses
      Lambda.df=Count.df
      Lambda.df[,"Adults"]=as.vector(t(Adults[,2:(tsteps.3+1)]/Adults[,1:tsteps.3]))
      Lambda.df[,"Pups"]=as.vector(t(Pups[,2:(tsteps.3+1)]/Pups[,1:tsteps.3]))
      for(idep in 1:2){
        Cur.df=Lambda.df
        if(idep==1){  #delete rows where lambda is NA in case of extirpation
          if(sum(is.na(Lambda.df[,"Adults"]))>0)Cur.df=Cur.df[-which(is.na(Lambda.df[,"Adults"])==1),]
          if(sum(Cur.df[,"Adults"]==Inf)>0)Cur.df=Cur.df[-which(Cur.df[,"Adults"]==Inf),]
        }
        else{
          if(sum(is.na(Lambda.df[,"Pups"]))>0)Cur.df=Cur.df[-which(is.na(Lambda.df[,"Pups"])==1),]          
          if(sum(Cur.df[,"Pups"]==Inf)>0)Cur.df=Cur.df[-which(Cur.df[,"Pups"]==Inf),]
        }
        for(iind in 1:4){
          dep=Dep.vars[idep]
          ind=Ind.vars[iind]
          fixed.formula=formula(eval(parse(text=paste(dep,'~',ind))))
          random.formula=formula(eval(parse(text=paste(dep,'~ Site'))))          
          cur.mod<-glmmLDTS(fixed.formula=fixed.formula,random.formula=random.formula,data=Cur.df,timecol="Time",group.vec="Site",distribution="Gaussian")          
          Result[isim,10+(idep-1)*4+iind]=sign(cur.mod$fixed.effects[2,"estimate"])*cur.mod$fixed.effects[2,"prob.t"]
          if(Result[isim,10+(idep-1)*4+iind]==0 & sign(cur.mod$fixed.effects[2,"estimate"])==-1)Result[isim,10+(idep-1)*4+iind]=-2
          
        }
      } 
      if((isim %% every)==0)write.csv(Result,file=cur.file)
    }
  }
}

#plot(sim.result$Svars$Ntot)
#plot(apply(sim.result$Fvars$F.fleet,2,'mean'))


