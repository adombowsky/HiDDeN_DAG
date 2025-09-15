power_simulate <- function(seed, n, k, epsilon) {
  sourceCpp("rcpp/mala.cpp")
  # seed = random number seed
  # n = sample size
  # k = number of categories in confounders
  # epsilon = MALA stepsize
  
  scores <- rep(0,4) # output vector
  
  # simulate confounders
  set.seed(seed)
  x1 <- sample(0:k,n,T)
  
  # simulate treatment
  ## conditional probabilities
  set.seed(seed)
  trt_prbs <- runif(min=.01,max=.99,n=k+1)
  ## data
  x2 <- rbinom(n=n,size=1,prob=trt_prbs[x1+1])
  
  # # simulate outcome
  x3 <- rep(0,n)
  set.seed(seed)
  out_prbs <- matrix(rbeta(shape1=2,shape2=15,n=2*(k+1)),nrow=k+1,ncol=2)
  #out_prbs <- matrix(runif(min=.01,max=.99,n=2*(k+1)),nrow=k+1,ncol=2)
  for (i in 1:n) {
    x3[i] <- rbinom(n=1,size=1,prob=out_prbs[x1[i]+1,x2[i]+1])
  }
  
  # create DF (for bnlearn)
  X_df <- data.frame(x_1 = factor(x1), x_2 = factor(x2), x_3 = factor(x3))
  
  # make two models
  g1 <- create_grids_and_table_child_rcpp(x3, matrix(x1))
  g2 <- create_grids_and_table_child_rcpp(x3, cbind(x1,x2))
  
  # lists
  ct <- list(g1$ct, g2$ct)
  mct <- list(g1$mct, g2$mct)
  par_grid <- list(g1$par_grid %>% as.matrix(), g2$par_grid %>% as.matrix())
  par_child_grid <- list(g1$par_child_grid %>% as.matrix(), g2$par_child_grid %>% as.matrix())
  
  # fit HiDDeN
  R <- 10000
  B <- 200
  set.seed(seed)
  pset_init <- sample(2,size=1,prob=c(.5,.5))
  fit <- mala_model_sampler_rcpp(R=R,
                                 k_j=1,
                                 ct = ct,
                                 mct = mct,
                                 par_grid = par_grid,
                                 par_child_grid = par_child_grid,
                                 b_j=1, rho_j=2,
                                 prior_probs = c(.5,.5),
                                 t_init = rep(1,2),
                                 pset_init = pset_init,
                                 step_size = epsilon,
                                 stops=R)
  colMeans(fit$t_accept[-(1:B),])
  scores[1] <- mean(fit$Pset[-(1:B)]==2)
  
  # comparison to score based techniques using bnlearn
  g_1 <- "[x_1][x_2|x_1][x_3|x_1]"
  g_1_dag <- model2network(g_1)
  g_2 <- "[x_1][x_2|x_1][x_3|x_1:x_2]"
  g_2_dag <- model2network(g_2)
  ## BIC
  scores[2] <- bnlearn::score(x=g_2_dag, data=X_df, type="bic") - bnlearn::score(x=g_1_dag, data=X_df, type="bic")
  ## AIC
  scores[3] <- bnlearn::score(x=g_2_dag, data=X_df, type="aic") - bnlearn::score(x=g_1_dag, data=X_df, type="aic")
  ## bde
  scores[4] <- bnlearn::score(x=g_2_dag, data=X_df, type="bde") - bnlearn::score(x=g_1_dag, data=X_df, type="bde")
  
  return(scores)
}