# preliminaries
library(foreach)
library(doParallel)
library(bnlearn)
library(Rcpp)
library(RcppArmadillo)
source("r/power_simulate.R")
#sourceCpp("rcpp/mala.cpp")
# parallel computing
n.cores <- 10
## create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
## register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# sample size
k <- c(4, 24, 99, 199) # number of categories in the counfounder
epsilon <- matrix(c(0.9,.3,
                    .65,.25,
                    .45,.15,
                    .4,.1), byrow=T, ncol=2)
S <- 100 # number of replications
# get simulations
sims <- foreach(h = 1:length(k)) %:%
  foreach(s=1:S,
          .packages = c("Rcpp","RcppArmadillo","bnlearn","dplyr")) %dopar% {
            power_simulate(seed=s,
                           n=200,
                           k=k[h],
                           epsilon = epsilon[h,])
          }
# post-process
# averages
avgs <- matrix(0, nrow = 4, ncol = 4)
for (j in 1:nrow(avgs)) {
  avgs[j,] <- colMeans(do.call(rbind, sims[[j]]))
}
sds <- matrix(0, nrow = 4, ncol = 4)
for (j in 1:nrow(avgs)) {
  sds[j,] <- apply(do.call(rbind, sims[[j]]), 2, sd)
}

# MAP estimator
HiDDeN <- matrix(0, nrow = S, ncol = 4)
for (j in 1:nrow(avgs)) {
  HiDDeN[,j] <- do.call(rbind, sims[[j]])[,1]
}
colMeans(HiDDeN>0.5)
# BIC estimator
bic_est <- matrix(0, nrow = S, ncol = 4)
for (j in 1:nrow(avgs)) {
  bic_est[,j] <- do.call(rbind, sims[[j]])[,2]
}
colMeans(bic_est>0)
# AIC estimator
aic_est <- matrix(0, nrow = S, ncol = 4)
for (j in 1:nrow(avgs)) {
  aic_est[,j] <- do.call(rbind, sims[[j]])[,3]
}
colMeans(aic_est>0)
# BDE estimator
BDE <- matrix(0, nrow = S, ncol = 4)
for (j in 1:nrow(avgs)) {
  BDE[,j] <- do.call(rbind, sims[[j]])[,4]
}
colMeans(BDE>0)

