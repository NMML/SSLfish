#plot selectivity, maturity, weight at age for SSL paper
source('c:/users/paul.conn/git/SSLfish/R/util_funcs.R')
a=c(1:10)
Sel=double.logistic(a,1.3,1,4,10)
Sel_SSL=double.logistic(a,1.3,1,2,8)
Mat=logistic(a,5,1)
Weight=logistic(a,5,0.7)

pdf("At_age_plots.pdf")

plot(a,Sel,type="l",lwd=3,lty=1,ylim=c(0,1),xlab="Age",ylab="Proportion or value",cex.axis=1.3,cex.lab=1.3)
lines(a,Sel_SSL,lwd=3,lty=2)
lines(a,Mat,lwd=3,lty=3)
lines(a,Weight,lwd=3,lty=5)
legend(4.8,.3,c("Fleet selectivity","SSL selectivity","Maturity","Weight"),lwd=c(3,3,3,3),lty=c(1,2,3,4),cex=1.3)

dev.off()