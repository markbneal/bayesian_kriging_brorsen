## Program written by B. Wade Brorsen, Oklahoma State University
## The program uses Monte Carlo generated data
## Bayesian Kriging is used to smooth spatial means across space
## The Stan part of the model runs about 25 seconds of System Time rather than 31 days
## Only the mean is smoothed across space
## presented at StanCon 2020

## Begin Stan portion
#now in separate .stan file

##Begin R program
library(mvtnorm)
library(parallel)
getOption("mc.cores", 1L) #Does nothing here since cores is specified below
library(rstan)
set.seed(97517)


nreps <- 40 #25 years is hard coded
nobs <- nreps*25
xx <- integer(length=nobs)
gg <- integer(length = nobs)
y <- vector("numeric",length = nobs)
error <- rnorm(nobs, mean=0,sd=10)
r_effect <- rnorm(25, mean=0,sd=10)
rho <- 6.75  ##has divergences with low values of rho


S <- matrix(runif(2*nreps),nreps,2)  # Create longitude-latitude coordinate
distance <- as.matrix(dist(S,upper=T,diag=T))  # Calculate distance
spatial_correlation = matrix(c(1:nreps),nreps,nreps)
sill <- 8; 

for (j in 1:nreps) {
  for (k in 1:nreps) {
    spatial_correlation[j,k]=sill*exp(-rho*distance[j,k]);
  }
}

#sill <- 8; #moved above

zz <- rnorm(nreps,mean=0,sd=1)
bb = rmvnorm(1,rep(46,nreps),spatial_correlation)
location =0;
for(n in 1:nobs) {
  II <- floor((n-1)/nreps)+1
  location = location + 1;
  if (location>nreps) {location=1}
  gg[n] <- II
  y[n] = bb[location]  + error[n]
}
replication <- rep(1:nreps,25)
xx <- cbind(xx)

mc_data <- list(N=NROW(y),M=25,L=nreps,y=y,dist=distance,g=gg, location=replication)

hfit <- stan(file="Bayesian_Kriging_Brorsen.stan", model_name="practice", data=mc_data, iter=2500, warmup=500, chains=4, cores=8,
             seed = 29473, control = list(adapt_delta = 0.80, max_treedepth = 15))


print(hfit, pars=c("a","sigmae","rho","sill"))
traceplot(hfit,pars=c("a","sigmae","rho","sill"),inc_warmup = TRUE)

list_of_draws <- extract(hfit)
print(names(list_of_draws))
pairs(hfit, pars=c("a","sigmae","rho","sill"))
