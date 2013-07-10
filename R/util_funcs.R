#inverse logit transform
invLogit <- function(x){
  1/(1+exp(-x))
}

# logit link function
logit <- function(x){
  if(any(x<0) | any(x>1)) stop("x must lie in (0,1)")
  log(x/(1-x))
}

#function to evaluate double logistic function for fishery selectivity
double.logistic<-function(A,eta1,eta2,kappa1,kappa2){
  s.a=1/(1+exp(-eta1*(A-kappa1)))*(1-1/(1+exp(-eta2*(A-kappa2))))
  s.a=s.a/max(s.a)
}

#function to compute logistic function for maturity, weight @ age
logistic<-function(A,a,b){
  m.a=(1+exp(-b*(A-a)))^{-1}
  m.a/max(m.a)
}

#Richards growth for SSL biomass
richards=function(age, A, m, S0, t){
  return((A^(1-m) - (A^(1-m) - S0^(1-m))*exp(-2*age*(1+m)/t))^(1/(1-m)))
}