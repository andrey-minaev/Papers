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
setwd("/Users/andrew/Dropbox/Research/Cost\ of\ privacy/Code/Estimations/Moraga-Gonzalez\ aka/Demand/Simmulated\ data")

load("data_simulated.RData")
dat <- data

invH0 <- function(s){
  #uniroot(function(r) dnorm(r,mean=0,sd=1)-r*(1-pnorm(r,mean=0,sd=1))-s, c(-10^200,10^200), tol = 1e-9, maxiter = 10000)$root
  if(s<100) uniroot(function(r) 0.57-r+incgam(exp(-r),0)-s, c(-100,100), tol = 1e-9)$root
  else -s
}

dat2 <- dat
dat <- NULL
for(i in unique(dat2$Consumer_ID)){
  dat[[i]] <- dat2[dat2$Consumer_ID == i,]
}



##################### Fixed coefficients model
#####################

LL <- function(theta, dataset){
  P <- NULL
  
  k=0
  for(i in unique(dataset$Consumer_ID)){
    delta <- t(theta[1:2]%*%t(dat[[i]][,c("Price","Quality")]))
    mu <- log(1+exp(theta[3]*dat[[i]]$Position))
    
    k=k+1
    P[k] <- exp((delta-mu)[dat[[i]]$booked==1])/sum(exp((delta-mu)))
  }
  print(theta)
  return(-sum(log(P)))
}



time <- Sys.time()
fit <- optim(rep(0,3), LL, dataset=dat2, method = c("L-BFGS-B"), #c("Nelder-Mead"),
             lower = c(rep(-Inf,2),0), upper = rep(Inf,3),
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


##################### Random coefficients model
#####################

d <- 1000 # number of random coefficients draws

coeff_unit <- mvrnorm(d,c(0,0),diag(c(1,1))) # random coefficients draws

LL <- function(theta, dataset){
  P <- NULL
  coeff <- matrix(NA,ncol=2,nrow=d)
  coeff[,1] <- theta[1] + coeff_unit[,1]*abs(theta[3])
  coeff[,2] <- theta[2] + coeff_unit[,2]*abs(theta[4])
  
  k=0
  for(i in unique(dataset$Consumer_ID)){
    delta <- matrix(NA,ncol=d,nrow=dim(dat[[i]])[1])
    delta <- dat[[i]]$Price%*%t(coeff[,1]) + dat[[i]]$Quality%*%t(coeff[,2])
    mu <- log(1+exp(theta[5]*dat[[i]]$Position))
    k=k+1
    P[k] <- mean(exp((delta-mu)[dat[[i]]$booked==1])/colSums(exp((delta-mu))))
  }
  print(theta)
  return(-sum(log(P)))
}



time <- Sys.time()
fit <- optim(c(0,0,1,1,0), LL, dataset=dat2, method = c("L-BFGS-B"), #c("Nelder-Mead"),
             #lower = c(rep(-Inf,2),rep(0,3)), upper = rep(Inf,5),
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
