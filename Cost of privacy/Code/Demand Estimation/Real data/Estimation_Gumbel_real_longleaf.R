library(MASS)
rm(list=ls())

dat2 <- read.csv("data_8192_ordered.csv")
dat2$const <- 1*(dat2$position>0)
dat <- NULL
for(i in unique(dat2$id)){
  dat[[i]] <- dat2[dat2$id == i,]
}

d <- 10000
coeff_unit <- mvrnorm(d,rep(0,6),diag(rep(1,6)))  

LL <- function(theta, dataset){
  P <- NULL
  coeff <- matrix(NA,ncol=6,nrow=d)
  coeff[,1] <- theta[1] + coeff_unit[,1]*abs(theta[7])
  coeff[,2] <- theta[2] + coeff_unit[,2]*abs(theta[8])
  coeff[,3] <- theta[3] + coeff_unit[,3]*abs(theta[9])
  coeff[,4] <- theta[4] + coeff_unit[,4]*abs(theta[10])
  coeff[,5] <- theta[5] + coeff_unit[,5]*abs(theta[11])
  coeff[,6] <- theta[6] + coeff_unit[,6]*abs(theta[12])
  
  k=0
  for(i in unique(dataset$id)){
    delta <- matrix(NA,ncol=d,nrow=dim(dat[[i]])[1])
    delta <- dat[[i]]$const%*%t(coeff[,1]) + dat[[i]]$price%*%t(coeff[,2]) + dat[[i]]$star%*%t(coeff[,3]) + dat[[i]]$review%*%t(coeff[,4]) + dat[[i]]$location%*%t(coeff[,5]) + dat[[i]]$brand%*%t(coeff[,6])
    mu <- log(1+exp(theta[13]*dat[[i]]$position))
    k=k+1
    P[k] <- mean(exp((delta-mu)[dat[[i]]$booked==1])/colSums(exp((delta-mu))))
  }
  iteration <<- iteration + 1
  if(iteration%%20==0){
    cat(paste("\n",iteration,"iterations. Time =",round(difftime(Sys.time(), time, units = "hours"),2),"hours","\n"))
  }
  if(iteration%%5==0){
    print(round(theta,3))
  }
  return(-sum(log(P)))
}

sink(file = "log_ordered_10000.txt", append = FALSE, type = c("output", "message"),
     split = FALSE)

time <- Sys.time()
iteration <- 0
fit <- optim(c(rep(0,6),rep(1,6),0), LL, dataset=dat2, method = c("L-BFGS-B"), #c("Nelder-Mead"),
             #lower = c(rep(-Inf,5),0), upper = rep(Inf,6),
             hessian=TRUE)#, control=list(trace=TRUE, REPORT=1,reltol=1e-3,abstol=1e-3))
Sys.time() - time

fisher_info<-solve(fit$hessian)
prop_sigma<-sqrt(diag(fisher_info))
upper<-fit$par+1.96*prop_sigma
lower<-fit$par-1.96*prop_sigma
interval<-data.frame(value=fit$par, upper=upper, lower=lower)
fit$par
is.positive.semi.definite(solve(-fit$hessian))
interval

sink()


