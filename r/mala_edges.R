library(Rcpp)
library(RcppArmadillo)
source("r/mala_within_gibbs.R")
sourceCpp("rcpp/mala.cpp")
create_grids_and_table_child <- function(x_j, x_Pa_j) {
  # x_j = vector of jth categorical variable
  # x_Pa_j = desired parent set for x_j, NONEMPTY
  
  # record number of categories
  k_j <- max(x_j)
  k_Pa_j <- apply(as.matrix(x_Pa_j), 2, max)
  
  # create new tables
  x_j_factored <- factor(x_j, levels = 0:k_j)
  x_Pa_j_factored <- Map(factor, as.list(as.data.frame(x_Pa_j)), lapply(k_Pa_j, function(x) 0:x))
  mct <- as.matrix(as.vector(table(x_Pa_j_factored)))
  ct <- as.matrix(as.vector(table(c(list(x_j_factored), x_Pa_j_factored))))
  par_grid <- expand.grid(lapply(k_Pa_j,function(x) 0:x))
  par_child_grid <- expand.grid(c(list(0:k_j), lapply(k_Pa_j,function(x) 0:x)))
  
  # returning
  return(list(ct=ct,
              mct=mct,
              par_child_grid=par_child_grid,
              par_grid=par_grid))
}

create_grids_and_table_root <- function(x_j) {
  # x_j = vector of jth categorical variable
  
  # record number of categories
  k_j <- max(x_j)
  k_Pa_j <- 0
  
  # create new tables
  mct <- matrix(length(x_j),nrow=1,ncol=1)
  ct <- as.matrix(as.vector(table(factor(x_j, levels = 0:k_j))))
  par_grid <- data.frame(Var1 = 0)
  par_child_grid <- expand.grid(list(0:k_j,0))
  
  
  # returning
  return(list(ct=ct,
              mct=mct,
              par_child_grid=par_child_grid,
              par_grid=par_grid))
}

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

update_edge <- function(t, z_minus, x_j, x_j_prime, X_minus, c_j, d_j) {
  # t = current value of latent prior means
  # z_minus = current edge indicator values for candidate parents (w/out j_prime)
  # x_j = variable for node j
  # x_j_prime = candidate parent j_prime
  # X_minus = data w/out x_j and x_j_prime
  
  M <- sum(z_minus) # size of current parent set w/out j_prime
  W <- sum(1-z_minus)
  
  if (M==0) {
    x_Pa_j_minus = NULL
  } else {
    x_Pa_j_minus = as.matrix(X_minus)[,z_minus==1]
  }
  
  # probability of no-edge
  if (M==0) {
    obj <- create_grids_and_table_root(x_j=x_j)
  } else {
    obj <- create_grids_and_table_child(x_j=x_j, x_Pa_j = x_Pa_j_minus)
  }
  
  w <- log( d_j + W  ) + log_node_score(t,
                                        obj$ct,
                                        obj$mct,
                                        obj$par_grid,
                                        obj$par_child_grid)
  
  # probability of an edge
  if (M==0) {
    obj <- create_grids_and_table_child(x_j=x_j, x_Pa_j=x_j_prime)
  } else {
    obj <- create_grids_and_table_child(x_j = x_j, x_Pa_j = cbind(x_Pa_j_minus, x_j_prime) )
  }
  
  m <- log(c_j + M) + log_node_score(t,
                                       obj$ct,
                                       obj$mct,
                                       obj$par_grid,
                                       obj$par_child_grid)
  
  # sample z_j_jprime
  
  z_j_j_prime <- sample(0:1, size = 1, prob = exp(c(w,m) - max(w,m)) )
  
  return(z_j_j_prime)
}

update_parents <- function(t, z, x_j, X_minus_j, c_j, d_j) {
  # update whole parent set
  for (j_prime in 1:length(z)) {
    z[j_prime] <- update_edge(t=t,
                              z_minus = z[-j_prime],
                              x_j = x_j,
                              x_j_prime = X_minus_j[,j_prime],
                              X_minus = X_minus_j[,-j_prime],
                              c_j=c_j,
                              d_j=d_j)
  }
  return(z)
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

mala_edge_sampler <- function(R,
                              x_j, X_minus,
                              b_j,rho_j,
                              c_j, d_j,
                              t_init=rep(1,max(x_j)+1),
                              step_size, stops=1) {
  # R = number of iterations
  # x_j = categorical variable for variable j
  # X_minus = all other categorical variables in Ca(j)
  # b, rho = hyperparameters
  # t_init = initial value of t
  # step_size = step size for MALA, assumed constant for all k (and graphs)
  
  
  ## initialize
  k_j <- max(x_j)
  z <- rbinom(ncol(as.matrix(X_minus)), 1, prob = (c_j/(c_j+d_j)))
  z_samps <- matrix(0, ncol = ncol(as.matrix(X_minus)), nrow = R)
  z_samps[1,] <- z
  t <- t_init
  t_samps <- matrix(0, nrow = R, ncol = k_j+1)
  t_samps[1,] <- t_init
  t_accept <- rep(1,k_j+1)
  t_accept_samps <- matrix(1, nrow = R, ncol = k_j+1)
  # u_samps <- matrix(0, nrow = R, ncol = nrow(par_grid))
  # u_samps[1,] <- u
  
  # configuring stops
  stops <- (1:(R/stops)) * stops
  
  ## run Gibbs
  for (r in 2:R) {
    # update latent variables
    if (sum(z)==0) {
      current_parents <- create_grids_and_table_root_rcpp(x_j=x_j)
    } else if (sum(z)==1) {
      current_parents <- create_grids_and_table_child_rcpp(x_j=x_j,
                                                      x_Pa_j=matrix(X_minus[,z==1],ncol=1))
    } else {
      current_parents <- create_grids_and_table_child_rcpp(x_j=x_j,
                                                           x_Pa_j=X_minus[,z==1])
    }
    mwg_update <- mala_within_gibbs_graph_rcpp(t=t,
                                          ct=current_parents$ct,
                                          mct=current_parents$mct,
                                          par_grid=current_parents$par_grid,
                                          par_child_grid=current_parents$par_child_grid,
                                          b=b_j,rho=rho_j,k=k_j,
                                          step_size=step_size,t_accept = rep(0,k_j+1))
    t <- mwg_update$t_update
    t_accept <- mwg_update$t_accept
    
    # update graph
    z <- update_parents_rcpp(t=t, z=z, x_j=x_j,
                        X_minus_j=X_minus,
                        c_j=c_j, d_j=d_j)
    
    ## exporting
    t_samps[r,] <- t
    #u_samps[r,] <- u
    t_accept_samps[r,] <- t_accept
    z_samps[r,] <- z
    
    # print stops
    if (r %in% stops){
      print(r)
    }
    
    ## debugging
    # print(t)
    # print(u)
  }
  return(list(t=t_samps,
              #u=u_samps,
              t_accept = t_accept_samps,
              z = z_samps))
}


