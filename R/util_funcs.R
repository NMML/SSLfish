#inverse logit transform
invLogit <- function(x){
  1/(1+exp(-x))
}

# logit link function
logit <- function(x){
  if(any(x<=0) | any(x>=1)) stop("x must lie in (0,1)")
  log(x/(1-x))
}

#function to evaluate double logistic function for fishery selectivity
double.logistic<-function(a,eta1,eta2,alpha1,alpha2){
  s.a=1/(1+exp(-eta1*(a-alpha1)))*(1-1/(1+exp(-eta2*(a-alpha2))))
  s.a=s.a/max(s.a)
}

#function to compute logistic function for maturity, weight @ age
logistic<-function(a,eta1,alpha1){
  m.a=(1+exp(-eta1*(a-alpha1)))^{-1}
  m.a/max(m.a)
}

BH.curve <- function(R0, h, SSB)

#compute unfished spawning biomass per recruit
unfished.bpr<-function(Maturity,Weight,M){
  n.age=length(Maturity)
  N=rep(1,n.age)
  for(i in 1:n.age)N[i]=exp(-M*i)
  N[n.age]=N[n.age]/(1-exp(-M))
  sum(N*Maturity*Weight)
}