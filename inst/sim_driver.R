#simdriver.R - Script to conduct SSL/fish/fisheries simulations
library(gtools)
nreps=1000
tsteps.1=100
tsteps.2=10
tsteps.3=35
tsteps.all=tsteps.1+tsteps.2+tsteps.3
n.yrs=tsteps.all+1
dem.opt=c("fec","surv")
effort.opt=c("rand","prop")
Result=data.frame(matrix(0,nreps*4,4))
colnames(Result=c("dem","effort","effect.pos","effect.neg"))
for(idem in 1:2){
  for(ieff in 1:2){
    dem=dem.opt[idem]
    eff=effort.opt[ieff]
    for(isim in 1:nreps){
      L.inits=initialize_pops()

      I
      Result[isim,1]=dem
      Result[isim,2]=eff
      Result[isim,3]=
      Result[isim,4]=
    }
  }
}

