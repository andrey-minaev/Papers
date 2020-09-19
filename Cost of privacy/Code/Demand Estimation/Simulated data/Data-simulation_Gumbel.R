library(MASS)
#while(!require(evd)){install.packages("evd")}
#while(!require(pracma)){install.packages("pracma")}
library(evd)
library(pracma)

#### Set working directory
setwd("/Users/andrew/Dropbox/Research/Cost\ of\ privacy/Code/Estimations/Moraga-Gonzalez\ aka/Demand/Simmulated\ data")

### Memory cleaning
rm(list=ls())

time <- Sys.time()
############### Main parameters of the model
n_hotels <- 200             # number of hotels
n_consumers <- 500        # number of consumers
hotels_per_consumer <- 30   # number of hotels shown to each consumer
k = 0.1                   # the marginal search cost per position



############### Generating hotels and consumers data

### Generating hotels distributional characteristics ~ Normal.

mu_hotels <- c(0,0)
mu_hotels[2] <- 4.877902
mu_hotels[1] <- 4.877902

sigma_p <- 1
sigma_q <- 1
corr_q_p <- 0.6

Sigma_hotels <- matrix(c(sigma_q**2, corr_q_p*sigma_q*sigma_p,
                         corr_q_p*sigma_q*sigma_p, sigma_p**2),
                       ncol = 2, nrow = 2)

### Generating hotels data

hotels <- mvrnorm(n = n_hotels, mu_hotels, Sigma_hotels, tol = 1e-6, empirical = FALSE)
hotels <- as.data.frame(hotels)
hotels <- cbind(hotels,seq(n_hotels))
names(hotels) <- c("Quality","Price", "ID")
hotels$Price <- exp(hotels$Price)/100
hotels$Quality <- exp(hotels$Quality)/100


### Generating consumers distributional characeristics. Alpha, beta ~ Normal

mu_consumers = c(-1,2)
sigma_alpha = .3
#sigma_alpha = 0
sigma_beta = .6
#sigma_beta = 0
corr_alpha_beta = 0
Sigma_consumers = matrix(c(sigma_alpha**2, corr_alpha_beta*sigma_alpha*sigma_beta,
                           corr_alpha_beta*sigma_alpha*sigma_beta, sigma_beta**2),
                         ncol = 2, nrow = 2)

### Generating consumers data data

consumers <- mvrnorm(n = n_consumers, mu_consumers, Sigma_consumers, tol = 1e-6, empirical = FALSE)
consumers <- as.data.frame(consumers)
names(consumers) <- c("alpha","beta")

#function invHO(s), inverse to H defined as function of s. Returns the reservation utility r to the search cost value s
invH0 <- function(s){
  #uniroot(function(r) dnorm(r,mean=0,sd=1)-r*(1-pnorm(r,mean=0,sd=1))-s, c(-10^200,10^200), tol = 1e-9, maxiter = 10000)$root
  if (s<100) uniroot(function(r) 0.57-r+incgam(exp(-r),0)-s, c(-100,100), tol = 1e-9)$root
  else -s
}

# Search cost cdf
Fc <- function(s,mu){
  (1-exp(-exp(-sapply(s, function(s) invH0(s))-mu)))/((1-exp(-exp(-sapply(s, function(s) invH0(s))))))
}

# Search cost generating function
search_cost_generator <- function(prob,mu){
  if (Fc(0,mu)<prob) uniroot(function(s) Fc(s,mu)-prob, c(0,100), tol = 1e-9)$root
  else 0
}


############### Combining consumers and hotels in queries dataframe.

### Generating matching. For each consumer random draw w/o replacement hotels_per_consumer hotels, which are displayed to this consumer.
match <- replicate(n_consumers,sample(seq(n_hotels), hotels_per_consumer, replace = FALSE, prob = NULL))

### Generating dataframe with positions and characteristics of drawn hotels for each consumer.
data <- data.frame(NULL)
for(i in 1:n_consumers){
  data <- rbind(data,c(i,rep(0,4)),cbind(rep(i,hotels_per_consumer),seq(hotels_per_consumer),hotels[match[,i],]))
}
names(data) <- c("Consumer_ID","Position","Quality","Price","Hotel_ID")


### Adding consumers' utilities u and w to the dataset.
# Generating Normal random shock
#epsilon <- rnorm((1+hotels_per_consumer)*n_consumers, 0, 1)
#epsilon0 <- rnorm(n_consumers, 0, 1)

# Generating TIEV random shock
epsilon <- rgumbel((1+hotels_per_consumer)*n_consumers, 0, 1)
epsilon0 <- rgumbel(n_consumers, 0, 1)


# Generating actual utility u, reservation utility r and w=min(u,r)
for(i in 1:dim(data)[1]){
  data$u0[i] <- epsilon0[data$Consumer_ID[i]]
  data$delta[i] <- sum(consumers[data$Consumer_ID[i],1:2]*data[i,c(4,3)])
  
  data$u[i] <- data$delta[i] + epsilon[i]
  data$r[i] <- data$delta[i] + invH0(search_cost_generator(runif(1),log(1+exp(k*data$Position[i]))))
  data$w[i] <- min(data$u[i],data$r[i])
}
# Outside option utilities correction
data$u[data$Position==0] <- data$u0[data$Position==0]
data$r[data$Position==0] <- data$u0[data$Position==0]
data$w[data$Position==0] <- data$u0[data$Position==0]

# Generating consumers purchasing decision. Each consumer purchases the hotel with highest reservaton calue w.

data$booked <- 0
for(i in 1:n_consumers){
  data[data$Consumer_ID==i,]$booked <- as.numeric(data[data$Consumer_ID==i,]$w==max(epsilon0[i],data[data$Consumer_ID==i,]$w))
}

# Generating consumers clicking desicion.
# data$clicked <- 0
# for(i in 1:n_consumers){
#   if(data$u0[data$Consumer_ID==i][1] < max(data$r[data$Consumer_ID==i])){
#     for(j in (((i-1)*(hotels_per_consumer+1))+1):(i*(hotels_per_consumer+1))){
#       if(data$r[j]>=max(data$w[data$Consumer_ID==i],data$u0[data$Consumer_ID==i][1])){
#         data$clicked[j] <- 1
#       }
#     }
#   }
# }  
# data$clicked[data$Position == 0] <- 1

save(data, file = "data_simulated.RData")

Sys.time() - time

#write.csv(data,file = "data_simulated.csv")








