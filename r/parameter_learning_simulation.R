### first study: binary variable messing with sparsity and number of categories ###

paramater_learning_simulation <- function(seed,n,k_par,R) {
  # dependencies
  library(Rcpp)
  library(RcppArmadillo)
  library(rstan)
  library(rstanarm)
  library(bnlearn)
  library(mgcv)
  source("r/mala_within_gibbs.R")
  # simulation study, one child and one parent but differing sample size/complexity
  ## seed = random number seed
  ## n = sample size
  ## k_par = number of categories in the parents (starting from 0)
  ## R = number of MCMC iterations (for Bayesian methods)
  
  rmse_vec <- c() # RMSE for HiDDeN and competitors
  
  # preliminaries
  ## simulating true conditional probabilities and data
  par_grid <- expand.grid(list(0:k_par))
  par_child_grid <- expand.grid(list(0:1,0:k_par))
  par_probs <- rep(1/(k_par+1),(k_par+1)) # category probabilities for the parent
  cond_probs <- matrix(0, nrow=k_par+1,ncol=2) # conditional probabilities for the child
  q <- 0
  for (h in 1:(k_par+1)) {
    q <- ifelse(h%%2==0, 1/3,2/3)
    cond_probs[h,] <- c(q,1-q)
  }
  set.seed(seed)
  x1 <- x2 <- rep(0,n)
  for (i in 1:n) {
    x1[i] <- sample(0:k_par,size=1,prob=par_probs)
  }
  for (i in 1:n) {
    x2[i] <- sample(0:1,size=1,prob=cond_probs[x1[i]+1,])
  }
  mct <- as.matrix(table(x1))
  ct <- as.matrix(as.vector(table(x2,x1)))
  dat <- data.frame(x_1=as.factor(x1),x_2=as.factor(x2))
  
  # run our MALA method
  step_size <- 5*10^(-1)*rep(1,2)
  HiDDeN_MALA <- mala_within_gibbs_sampler(R=R,
                                           ct=ct,mct=mct,
                                           par_grid=par_grid,par_child_grid = par_child_grid,
                                           b=1,rho=2,k=1,
                                           step_size=step_size)
  HiDDeN_cond_probs_samps <- classifier_MCMC_samples(t_samps=HiDDeN_MALA$t, 
                                                     ct=ct, mct=mct, 
                                                     par_child_grid = par_child_grid, par_grid = par_grid)
  HiDDeN_posterior_mean <- colMeans(HiDDeN_cond_probs_samps)
  HiDDeN_posterior_mean_success <- HiDDeN_posterior_mean[seq(2,length(HiDDeN_posterior_mean),by=2)] 
  rmse_vec[1] <- sqrt( sum((HiDDeN_posterior_mean_success-cond_probs[,2])^2)/(k_par+1) ) # HiDDeN RMSE
  
  # competitor: logistic regression (Fisher-Scoring)
  ct_df <- as.data.frame(matrix(ct, nrow=k_par+1, byrow=T))
  ct_df$parent <- as.factor(1:nrow(ct_df))
  log_reg <- glm(cbind(V2,V1)~parent,
                 data=ct_df,
                 family=binomial(link="logit"))
  rmse_vec[2] <- sqrt( sum((log_reg$fitted.values-cond_probs[,2])^2)/(k_par+1) ) # logistic regression RMSE
  
  # competitor: Bayesian logistic regression
  log_reg_bayes <- stan_glm(cbind(V2,V1)~parent,
                            data=ct_df,
                            family=binomial(link="logit"),
                            chains=4,
                            iter=R,
                            seed=seed,
                            refresh=0 # no printing iterations
  )
  rmse_vec[3] <- sqrt( sum((log_reg_bayes$fitted.values-cond_probs[,2])^2)/(k_par+1) ) # Bayesian logistic regression RMSE
  
  # competitor: GAM
  log_reg_gam <- gam(cbind(V2,V1)~parent,
                     data=ct_df,
                     family=binomial(link="logit"))
  rmse_vec[4] <- sqrt( sum((log_reg_gam$fitted.values-cond_probs[,2])^2)/(k_par+1) ) # GAM RMSE
  
  # competitor: MLEs
  MLEs <- ct_df$V2/mct
  rmse_vec[5] <- sqrt( sum((MLEs-cond_probs[,2])^2)/(k_par+1) ) # MLEs RMSE
  
  # competitor: bnlearn
  ## creating simple chain structure
  chain_net <- model2network("[x_1][x_2|x_1]")
  ## marginal MLEs
  cpt_x1 <- matrix(mct/n, ncol = k_par+1,
                   dimnames = list(NULL, as.factor(1:(k_par+1))))
  
  
  cpt_x2 <- array(as.vector(t(as.matrix(ct_df[,1:2]/(rowSums(ct_df[,1:2]))))), 
                  dim = c(2,k_par+1), 
                  dimnames = list(
                    "x_2" = c("1", "2"),
                    "x_1" = as.factor(1:(k_par+1))
                  ))
  bde <- bn.fit(chain_net, dat, method = "bayes")
  rmse_vec[6] <- sqrt( sum((bde$x_2$prob[2,]-cond_probs[,2])^2)/(k_par+1) ) # bde RMSE
  
  return(rmse_vec)
}



