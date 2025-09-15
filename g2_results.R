# preliminaries
library(foreach)
library(doParallel)
source("r/g2_simulations.R")
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
N <- c(50, 75, 100, 150)
S <- 100 # number of replications
# get simulations
sims <- foreach(i = 1:length(N)) %:%
  foreach(s=1:S,
          .packages = c("bnlearn")) %dopar% {
            g2_simulations(seed=s,
                           n=N[i])
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
