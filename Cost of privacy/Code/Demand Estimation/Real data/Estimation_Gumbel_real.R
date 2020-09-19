#while(!require(Rfast)){install.packages("Rfast")}
#while(!require(foreach)){install.packages("foreach")}
#while(!require(doParallel)){install.packages("doParallel")}
#while(!require(matrixcalc)){install.packages("matrixcalc")}
#while(!require(optimx)){install.packages("optimx")}
#library(Rfast)
#library(foreach)
#library(doParallel)
library(matrixcalc)
library(optimx)

# cores=detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# clusterExport(cl, c("dnorm", "pnorm", "qnorm", "rnorm"))
# registerDoParallel(cl)

rm(list=ls())
setwd("/Users/andrew/Dropbox/Research/Cost\ of\ privacy/Code/Estimations/Moraga-Gonzalez\ aka/Demand/Real\ data/Fixed-coeff")
dat <- read.csv("data_8192_ordered.csv")
dat <- dat[,c("id","star","review","brand","location","position","price","clicked","booked")]

invH0 <- function(s){
  #uniroot(function(r) dnorm(r,mean=0,sd=1)-r*(1-pnorm(r,mean=0,sd=1))-s, c(-10^200,10^200), tol = 1e-9, maxiter = 10000)$root
  if(s<100) uniroot(function(r) 0.57-r+incgam(exp(-r),0)-s, c(-100,100), tol = 1e-9)$root
  else -s
}

dat2 <- dat
dat <- NULL
for(i in unique(dat2$id)){
  dat[[i]] <- dat2[dat2$id == i,]
}



##################### Fixed coefficients NO constant, WITH search
#####################

LL <- function(theta, dataset){
  P <- NULL
  #inv <- sapply(log(1+exp(theta[11]+theta[12]*(0:38))),function(x) invH0(x))
  
  k=0
  for(i in unique(dataset$id)){
    delta <- t(theta[1:5]%*%t(dat[[i]][,c("price","star","review","brand","location")]))
    mu <- log(1+exp(theta[6]*dat[[i]]$position))
    
    k=k+1
    P[k] <- exp((delta-mu)[dat[[i]]$booked==1])/sum(exp((delta-mu)))
  }
  print(theta)
  return(-sum(log(P)))
}



time <- Sys.time()
fit <- optim(rep(0,6), LL, dataset=dat2, method = c("L-BFGS-B"), #c("Nelder-Mead"),
             #lower = c(rep(-Inf,5),0), upper = rep(Inf,6),
             hessian=TRUE, control=list(trace=TRUE, REPORT=1,reltol=1e-3,abstol=1e-3))
Sys.time() - time

fisher_info<-solve(fit$hessian)
prop_sigma<-sqrt(diag(fisher_info))
upper<-fit$par+1.96*prop_sigma
lower<-fit$par-1.96*prop_sigma
interval<-data.frame(value=fit$par, upper=upper, lower=lower)
fit$par
is.positive.semi.definite(solve(-fit$hessian))
interval


##################### Fixed coefficients with constant, WITH search
#####################


LL <- function(theta, dataset){
  P <- NULL
  
  k=0
  for(i in unique(dataset$id)){
    delta <- t(theta[1:6]%*%t(cbind(rep(1,dim(dat[[i]])[1]),dat[[i]][,c("price","star","review","brand","location")])))
    mu <- log(1+exp(theta[7]*dat[[i]]$position))
    
    k=k+1
    P[k] <- exp((delta-mu)[dat[[i]]$booked==1])/sum(exp((delta-mu)))
  }
  print(theta)
  return(-sum(log(P)))
}



time <- Sys.time()
fit <- optim(rep(0,7), LL, dataset=dat2, method = c("L-BFGS-B"), #c("Nelder-Mead"),
             #lower = c(rep(-Inf,5),0), upper = rep(Inf,6),
             hessian=TRUE, control=list(trace=TRUE, REPORT=1,reltol=1e-3,abstol=1e-3))
Sys.time() - time

fisher_info<-solve(fit$hessian)
prop_sigma<-sqrt(diag(fisher_info))
upper<-fit$par+1.96*prop_sigma
lower<-fit$par-1.96*prop_sigma
interval<-data.frame(value=fit$par, upper=upper, lower=lower)
fit$par
is.positive.semi.definite(solve(-fit$hessian))
interval


##################### Random coefficients, WITH search
#####################

d <- 1000
#d <- 5

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
    delta <- dat[[i]]$price%*%t(coeff[,1]) + dat[[i]]$star%*%t(coeff[,2]) + dat[[i]]$review%*%t(coeff[,3]) + dat[[i]]$location%*%t(coeff[,4]) + dat[[i]]$brand%*%t(coeff[,5])
    mu <- log(1+exp(theta[13]*dat[[i]]$position))
    k=k+1
    P[k] <- mean(exp((delta-mu)[dat[[i]]$booked==1])/colSums(exp((delta-mu))))
  }
  print(theta)
  return(-sum(log(P)))
}


time <- Sys.time()
fit <- optim(c(rep(0,6),rep(1,6),0), LL, dataset=dat2, method = c("L-BFGS-B"), #c("Nelder-Mead"),
             #lower = c(rep(-Inf,5),0), upper = rep(Inf,6),
             hessian=TRUE, control=list(trace=TRUE, REPORT=1,reltol=1e-3,abstol=1e-3))
Sys.time() - time

fisher_info<-solve(fit$hessian)
prop_sigma<-sqrt(diag(fisher_info))
upper<-fit$par+1.96*prop_sigma
lower<-fit$par-1.96*prop_sigma
interval<-data.frame(value=fit$par, upper=upper, lower=lower)
fit$par
is.positive.semi.definite(solve(-fit$hessian))
interval


