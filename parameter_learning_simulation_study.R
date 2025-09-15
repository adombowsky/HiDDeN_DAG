# preliminaries
library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doParallel)
# r functions
source("r/parameter_learning_simulation.R")
# parallel computing
n.cores <- 10
## create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
## register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# in parallel
N <- 100 # sample size
K <- c(1,2,4,9)
S <- 50
sims <- foreach(i = 1:length(N)) %:%
  foreach(k = 1:length(K)) %:%
  foreach(s=1:S,
          .packages = c("Rcpp", "RcppArmadillo",
                        "rstan","rstanarm",
                        "bnlearn","mgcv")) %dopar% {
          paramater_learning_simulation(seed=s,
                                        n=N[i],
                                        k_par=K[k],
                                        R=10000)
            
                        }

# returning column means
num_N <- length(N)
num_k_par <- length(K)

# configuring averages
num_metrics <- length(sims[[1]][[1]][[1]])
metric_names <- paste0("Metric_", 1:num_metrics)
averaged_results_by_k_par <- vector("list", num_k_par)
names(averaged_results_by_k_par) <- paste0("k_par_val_", K)

for (k_idx in 1:num_k_par) {
  # rows correspond to N values, columns to metrics
  avg_matrix_for_k <- matrix(NA, nrow = num_N, ncol = num_metrics)
  colnames(avg_matrix_for_k) <- metric_names # Apply metric names to columns
  rownames(avg_matrix_for_k) <- paste0("N_val_", N) # Apply N values to rows
  
  for (i_idx in 1:num_N) {
    results_for_ik_list <- sims[[i_idx]][[k_idx]]
    combined_results_matrix <- do.call(rbind, results_for_ik_list)
    avg_metrics_for_ik <- round(colMeans(combined_results_matrix),3)
    avg_matrix_for_k[i_idx, ] <- avg_metrics_for_ik
  }
  averaged_results_by_k_par[[k_idx]] <- avg_matrix_for_k
}

# configuring sds
std_dev_results_by_k_par <- vector("list", num_k_par)
names(std_dev_results_by_k_par) <- paste0("k_par_val_", K)

# Populate each matrix in the list
for (k_idx in 1:num_k_par) {
  sd_matrix_for_k <- matrix(NA, nrow = num_N, ncol = num_metrics)
  colnames(sd_matrix_for_k) <- metric_names # Apply metric names to columns
  rownames(sd_matrix_for_k) <- paste0("N_val_", N) # Apply N values to rows
  
  for (i_idx in 1:num_N) {
    results_for_ik_list <- sims[[i_idx]][[k_idx]]
    combined_results_matrix <- do.call(rbind, results_for_ik_list)
    std_dev_metrics_for_ik <- round(apply(combined_results_matrix, 2, sd),3)
    sd_matrix_for_k[i_idx, ] <- std_dev_metrics_for_ik
  }
  std_dev_results_by_k_par[[k_idx]] <- sd_matrix_for_k
}