g1_simulations <- function(seed,
                           n,
                           k_1=4,
                           k_2=2,
                           k_3=1,
                           R=20000,
                           B=200,
                           step_size=1.5*10^(-1)*rep(1,2)) {
  # dependencies
  source("r/mala_within_gibbs.R")
  source("r/mala_graph.R")
  library(bnlearn)
  
  score_vec <- c() # saving scores
  
  # preliminaries
  ## simulating true conditional probabilities and data
  par_probs <- rep(1/(k_2+1),(k_2+1)) # category probabilities for x_2
  root_probs <- rep(1/(k_1+1),k_1+1) # category probabilities for x_1
  cond_probs <- matrix(0, nrow = (k_1+1), ncol = (k_2+1)) # conditional probabilities for the child
  set.seed(seed)
  for (h in 1:(k_1+1)) {
    for (i in 1:(k_2+1)) {
      cond_probs[h,i] <- runif(n=1,min=0,max=1)
    }
  }
  x1 <- x2 <- x3 <- rep(0,n)
  for (i in 1:n) {
    x1[i] <- sample(0:k_1,size=1,prob=root_probs)
  }
  for (i in 1:n) {
    x2[i] <- sample(0:k_2,size=1,prob=par_probs)
  }
  for (i in 1:n) {
    x3[i] <- rbinom(n=1,size=1,p=1-cond_probs[x1[i]+1, x2[i]+1])
  }
  dat <- data.frame(x_1=as.factor(x1),
                    x_2=as.factor(x2), 
                    x_3=as.factor(x3))
  # under graph 1
  mct_g1 <- as.matrix(as.vector(table(x1,x2)))
  ct_g1 <- as.matrix(as.vector(table(x3,x1,x2)))
  par_grid_g1 <- expand.grid(list(0:k_1,0:k_2))
  par_child_grid_g1 <- expand.grid(list(0:1,0:k_1,0:k_2))
  
  # under graph 2
  mct_g2 <- as.matrix(table(x2))
  ct_g2 <- as.matrix(as.vector(table(x3,x2)))
  par_grid_g2 <- expand.grid(list(0:k_2))
  par_child_grid_g2 <- expand.grid(list(0:1,0:k_2))
  
  # fit out model
  fit <- mala_graph_sampler(R=R,
                            ct=list(ct_g1,ct_g2),
                            mct = list(mct_g1, mct_g2),
                            par_grid = list(par_grid_g1, par_grid_g2),
                            par_child_grid = list(par_child_grid_g1, par_child_grid_g2),
                            b=1,rho=2,k=1,
                            prior_probs = c(1/2,1/2),
                            step_size = step_size)
  score_vec[1] <- p_g <- mean(fit$pset[-(1:B)]==1)
  
  # comparison to score based techniques using bnlearn
  g_1 <- "[x_1][x_2][x_3|x_1:x_2]"
  g_1_dag <- model2network(g_1)
  g_2 <- "[x_1][x_2][x_3|x_2]"
  g_2_dag <- model2network(g_2)
  ## BIC
  score_vec[2] <- bnlearn::score(x=g_1_dag, data=dat, type="bic") - bnlearn::score(x=g_2_dag, data=dat, type="bic")
  
  ## AIC
  score_vec[3] <- bnlearn::score(x=g_1_dag, data=dat, type="aic") -bnlearn::score(x=g_2_dag, data=dat, type="aic")
  
  ## bde
  score_vec[4] <- bnlearn::score(x=g_1_dag, data=dat, type="bde") -bnlearn::score(x=g_2_dag, data=dat, type="bde")
  return(score_vec)
}

