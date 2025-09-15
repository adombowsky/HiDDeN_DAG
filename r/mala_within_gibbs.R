### implements the MALA-within-Gibbs sampler for an arbitrary network ###
source("r/mala.R")


## a single Gibbs iteration ##
mala_within_gibbs <- function(t,u,
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
  
  ## first, MALA on all categories of the child
  for (x in 1:(k+1)) {
    t_x_mala <- mala_update(t_current = t[x], step_size = step_size[x],
                            b=b, rho=rho, k=k, 
                            N=ct[par_child_grid[,1]==(x-1)],
                            U=u)
    t[x] <- t_x_mala$next_state
    t_accept[x] <- t_x_mala$acceptance
  }
  
  ## next, Uniform sampling on all categories of the parents
  for (x_par in 1:nrow(par_grid)) {
    u[x_par] <- beta_update(sum(t), mct[x_par])
  }
  return(list(t_update = t, u_update=u, t_accept = t_accept))
}

## a Gibbs sampler for any node ##
mala_within_gibbs_sampler <- function(R,
                                      ct,mct,par_grid,par_child_grid,
                                      b,rho,k,
                                      step_size) {
  # R = number of iterations
  # ct = parent-child cell counts
  # mct = parent cell counts (marginal cell counts)
  # par_grid = parent grid
  # par_child_grid = parent-child grid
  # b, rho = hyperparameters
  # k = child categories
  # step_size = step size for MALA, assumed constant for all k
  
  ## initialize
  t <- rexp(k+1)
  t_samps <- matrix(0, nrow = R, ncol = k+1)
  t_samps[1,] <- t
  t_accept <- rep(1,k+1)
  t_accept_samps <- matrix(1, nrow = R, ncol = k+1)
  u <- runif(nrow(par_grid))
  u_samps <- matrix(0, nrow = R, ncol = nrow(par_grid))
  u_samps[1,] <- u
  
  ## run Gibbs
  for (r in 2:R) {
    mwg_update <- mala_within_gibbs(t,u,
                                    ct,mct,par_grid,par_child_grid,
                                    b,rho,k,
                                    step_size,t_accept = rep(0,k+1))
    t <- mwg_update$t_update
    u <- mwg_update$u_update
    t_accept <- mwg_update$t_accept
    
    ## exporting
    t_samps[r,] <- t
    u_samps[r,] <- u
    t_accept_samps[r,] <- t_accept
    
    ## debugging
    # print(t)
    # print(u)
  }
  return(list(t=t_samps, u=u_samps,t_accept = t_accept_samps))
}

# ### example
# R <- 10^4
# step_size = 0.2*rep(1,k+1)
# fit <- mala_within_gibbs_sampler(R=R,
#                                  ct=ct, mct=mct,
#                                  par_grid = par_grid,
#                                  par_child_grid=par_child_grid,
#                                  b=b,rho=rho,k=k,
#                                  step_size=step_size)
# colMeans(fit$t_accept) # acceptance probability
# plot(fit$t[,1])
# plot(fit$t[,2])
# coda::effectiveSize(fit$t)

## convert to conditional probabilities

classifier_MCMC <- function(t, ct, mct, par_child_grid, par_grid) {
  # t = a single MCMC sample of t
  # ct = parent cell counts (marginal cell counts)
  # par_child_grid = parent-child grid
  # par_grid = parent_grid
  classifier <- rep(0, nrow(par_child_grid))
  beta_j <- sum(t)
  for (h in 1:nrow(par_child_grid)) {
    classifier[h] <- (t[par_child_grid[h,1]+1] + ct[h])/(beta_j + 
    mct[apply(par_grid,1,function(x) all(x==par_child_grid[h,-1]))])
  }
  return(classifier)
}

classifier_MCMC_samples <- function(t_samps, ct, mct, par_child_grid, par_grid) {
  classifier_samps <- t(apply(t_samps,1,classifier_MCMC,
                            ct=ct,mct=mct,par_child_grid=par_child_grid,par_grid=par_grid))
}
