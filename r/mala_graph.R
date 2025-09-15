source("r/mala_within_gibbs.R")

log_node_score <- function(t,ct,mct,par_grid,par_child_grid) {
  beta_j <- sum(t)
  A <- 0 # normalizing constant
  for (x in 1:nrow(par_grid)) {
    A <- A + lgamma(beta_j) - lgamma(mct[x] + beta_j)
  }
  B <- 0 # node-parent kernel
  for (x in 1:nrow(par_child_grid)) {
    B <- B + lgamma( ct[x] + t[par_child_grid[x,1]+1] ) - lgamma(t[par_child_grid[x,1]+1])
  }
  return(A+B)
}

parent_update <- function(t, ct, mct, par_grid, par_child_grid,
                          prior_probs = rep(1/length(ct), length(ct))) {
  # NOTE: for non-root hypotheses
  # t = vector of length k_j
  # ct = list of parent-child cell counts for each DAG
  # mct = list of parent cell counts for each DAG
  # par_grid = list of parent grids for each DAG
  # par_child_grid = list of parent-child grids for each DAG
  # prior_probs = vector of prior probabilities for each DAG
  
  M <- length(ct) # possible DAGs
  l_probs <- rep(0,M)
  for (m in 1:M) {
    l_probs[m] <- log(prior_probs[m]) +
      log_node_score(t,ct[[m]],mct[[m]],par_grid[[m]],par_child_grid[[m]])
  }
  parent_set <- sample(1:M, size=1, replace=T, prob = exp(l_probs - max(l_probs)))
  return(parent_set)
}

mala_within_gibbs_graph <- function(t,
                              ct,mct,par_grid,par_child_grid,
                              b,rho,k,
                              step_size,t_accept = rep(0,k+1)){
  # ct = parent-child cell counts
  # mct = parent cell counts (marginal cell counts)
  # par_grid = parent grid
  # par_child_grid = parent-child grid
  # b, rho = hyperparameters
  # k = child categories
  # step_size = step size for MALA, a vector of length k+1
  # t_accept = indicators for t_acceptance
  
  ## first, Uniform sampling on all categories of the parents
  u <- rep(0, nrow(par_grid))
  for (x_par in 1:nrow(par_grid)) {
    u[x_par] <- beta_update(sum(t), mct[x_par])
  }
  
  ## next, MALA on all categories of the child
  for (x in 1:(k+1)) {
    t_x_mala <- mala_update(t_current = t[x], step_size = step_size[x],
                            b=b, rho=rho, k=k, 
                            N=ct[par_child_grid[,1]==(x-1)],
                            U=u)
    t[x] <- t_x_mala$next_state
    t_accept[x] <- t_x_mala$acceptance
  }
  return(list(t_update = t, u_update=u, t_accept = t_accept))
}

mala_graph_sampler <- function(R,
                               ct,mct,par_grid,par_child_grid,
                               b,rho,k,
                               prior_probs,
                               step_size) {
  # R = number of iterations
  # ct = list of parent-child cell counts
  # mct = list of parent cell counts (marginal cell counts)
  # par_grid = list of parent grids
  # par_child_grid = list of parent-child grids
  # b, rho = hyperparameters
  # k = child categories
  # prior_probs = prior probabilties of graphs
  # step_size = step size for MALA, assumed constant for all k (and graphs)
  
  
  ## initialize
  pset <- sample(1:length(ct), size = 1, replace = T, prob = prior_probs)
  pset_samps <- rep(pset,R)
  t <- rexp(k+1)
  t_samps <- matrix(0, nrow = R, ncol = k+1)
  t_samps[1,] <- t
  t_accept <- rep(1,k+1)
  t_accept_samps <- matrix(1, nrow = R, ncol = k+1)
  # u_samps <- matrix(0, nrow = R, ncol = nrow(par_grid))
  # u_samps[1,] <- u
  
  ## run Gibbs
  for (r in 2:R) {
    # update latent variables
    mwg_update <- mala_within_gibbs_graph(t=t,
                                    ct=ct[[pset]],mct=mct[[pset]],par_grid=par_grid[[pset]],par_child_grid=par_child_grid[[pset]],
                                    b=b,rho=rho,k=k,
                                    step_size=step_size,t_accept = rep(0,k+1))
    t <- mwg_update$t_update
    t_accept <- mwg_update$t_accept
    
    # update graph
    pset <- parent_update(t=t, ct=ct, mct=mct, par_grid=par_grid, par_child_grid=par_child_grid,
                          prior_probs = prior_probs)
    
    ## exporting
    t_samps[r,] <- t
    #u_samps[r,] <- u
    t_accept_samps[r,] <- t_accept
    pset_samps[r] <- pset
    
    ## debugging
    # print(t)
    # print(u)
  }
  return(list(t=t_samps, 
              #u=u_samps,
              t_accept = t_accept_samps, pset = pset_samps))
}