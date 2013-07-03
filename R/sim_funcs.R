
###
### Functions to increment populations by one year
###

### Sea lion functions
library(compiler)

inc.N.ssl <- function(A.list, N.ssl.list){
  
}

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