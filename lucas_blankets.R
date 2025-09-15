source("r/mala_edges.R")
library(dplyr)
library(bnlearn)
library(glmnet)
library(ggplot2)
library(cowplot)
library(latex2exp)

## part 0: functions needed
## compute the SHD for HC and TABU
bnl_shd <- function(lasso_parents, z_true) {
  return(sum(lasso_parents!=z_true))
}
## obtain coefficient names for lasso
check_coefs <- function(vec, prefixes) {
  var_names <- unique(substr(vec, 1, nchar(vec)-1))
  return(as.numeric(prefixes %in% var_names))
}

# shd matrix
shds <- matrix(0,nrow=6,ncol=4)

## LUCAS0
dat <- read.csv("data/lucas.csv")
dat_df <- as.data.frame(dat) %>%
  mutate(across(.fns = as.factor)) ## for lasso
colnames(dat_df) <- paste0("V", 1:ncol(dat))


### Smoking
p_outcome <- which(colnames(dat)=="Smoking")
z_true <- as.numeric( colnames(dat %>% select(-Smoking)) %in% 
                        c("Anxiety","Peer_Pressure", "Yellow_Fingers", "Lung_cancer", "Genetics")  )
## Gibbs
R <- 10000
B <- 200 # burn-in
epsilon <- c(0.3,0.3)
set.seed(1)
fit <- mala_edge_sampler_rcpp(R=R,
                              x_j=dat$Smoking,
                              X_minus = dat %>% select(-Smoking) %>% as.matrix(),
                              b_j=1, rho_j=2,
                              c_j=1, d_j=1,
                              step_size=epsilon,
                              t_init = rep(1,2),
                              stops=200)
colMeans(fit$t_accept[-(1:B),])
saveRDS(fit, "lucas_results/smoking.rds")
### extract parent samples
parent_strings <- apply(fit$z[-(1:B),],1,paste,collapse="")
table(parent_strings)
### MAP estimate
MAP <- as.numeric(strsplit(names(which.max(table(parent_strings))),split="")[[1]])
### posterior Hamming distance
#shd_samples <- apply(fit$z[-(1:B),],1,function(y) sum(y!=z_true))
shds[1,1] <- bnl_shd(MAP,z_true)
### log posterior
log_posterior_smoking <- log_posterior_samples(t=fit$t[-(1:B),],
                                              z=fit$z[-(1:B),],
                                              x_j = dat$Smoking,
                                              X_minus = dat %>% select(-Smoking) %>% as.matrix(),
                                              b_j=1, rho_j=2,
                                              c_j=1, d_j=1)
lp_df <- data.frame(iter = (B+1):R,
                    smoking = log_posterior_smoking)
trace_smoking<- ggplot(lp_df, aes(x=iter,y=smoking)) + geom_line() + 
  xlab("Iteration") + ylab("Log-Posterior") + labs(title = "Smoking Log-Posterior MCMC Samples")  + theme_bw()
### lasso
formula_lasso <- as.formula(paste0("V", p_outcome, " ~ ."))
mm <- model.matrix(formula_lasso, data = dat_df)[,-1] 
cv_lasso <- cv.glmnet(x=mm, 
                      y=dat_df$V1, 
                      family = "multinomial",
                      type.multinomial="grouped",
                      alpha = 1) ## cross validation
coef(cv_lasso, s = "lambda.1se")
nz_ind <- predict(cv_lasso, s = "lambda.1se", type = "nonzero") # get nonzero coefficients
if (nz_ind$lambda.1se[1]==1) {
  nz_ind <- nz_ind %>% slice(-1) # drop intercept
}
coef_names <- rownames(coef(cv_lasso,s="lambda.1se")[[1]])
nz_coef_names <- coef_names[nz_ind$lambda.1se]
lasso_parents <- check_coefs(vec = nz_coef_names, prefixes=colnames(dat_df[,-p_outcome]))
shds[1,2] <- bnl_shd(lasso_parents=lasso_parents, z_true = z_true)
### bnlearn
gs_fit <- learn.mb(dat_df, paste0("V",toString(p_outcome)), method="gs")
gs_parents <- as.numeric(colnames(dat_df[,-p_outcome])%in% gs_fit)
shds[1,3] <- bnl_shd(gs_parents, z_true)
iamb_fit <- learn.mb(dat_df, paste0("V",toString(p_outcome)), method="iamb")
iamb_parents <- as.numeric(colnames(dat_df[,-p_outcome])%in% iamb_fit)
shds[1,4] <- bnl_shd(iamb_parents, z_true)



### Genetics
p_outcome <- which(colnames(dat)=="Genetics")
z_true <- as.numeric( colnames(dat %>% select(-Genetics)) %in% 
                        c("Lung_cancer","Attention_Disorder", "Smoking")  )
## Gibbs
R <- 10000
B <- 200 # burn-in
epsilon <- c(0.47,0.18)
set.seed(1)
fit <- mala_edge_sampler_rcpp(R=R,
                              x_j=dat$Genetics,
                              X_minus = dat %>% select(-Genetics) %>% as.matrix(),
                              b_j=1, rho_j=2,
                              c_j=1, d_j=1,
                              step_size=epsilon,
                              t_init = rep(1,2),
                              stops=200)
colMeans(fit$t_accept[-(1:B),])
saveRDS(fit, "lucas_results/genetics.rds")
### extract parent samples
parent_strings <- apply(fit$z[-(1:B),],1,paste,collapse="")
table(parent_strings)
### MAP estimate
MAP <- as.numeric(strsplit(names(which.max(table(parent_strings))),split="")[[1]])
### posterior Hamming distance
#shd_samples <- apply(fit$z[-(1:B),],1,function(y) sum(y!=z_true))
shds[2,1] <- bnl_shd(MAP,z_true)
## log-posterior
log_posterior_genetics <- log_posterior_samples(t=fit$t[-(1:B),],
                                               z=fit$z[-(1:B),],
                                               x_j = dat$Genetics,
                                               X_minus = dat %>% select(-Genetics) %>% as.matrix(),
                                               b_j=1, rho_j=2,
                                               c_j=1, d_j=1)
lp_df <- data.frame(iter = (B+1):R,
                    genetics = log_posterior_genetics)
trace_genetics <- ggplot(lp_df, aes(x=iter,y=genetics)) + geom_line() + 
  xlab("Iteration") + ylab("Log-Posterior") + labs(title = "Genetics Log-Posterior MCMC Samples")  + theme_bw()
### lasso
formula_lasso <- as.formula(paste0("V", p_outcome, " ~ ."))
mm <- model.matrix(formula_lasso, data = dat_df)[,-1] 
cv_lasso <- cv.glmnet(x=mm, 
                      y=dat_df$V5, 
                      family = "multinomial",
                      type.multinomial="grouped",
                      alpha = 1) ## cross validation
coef(cv_lasso, s = "lambda.1se")
nz_ind <- predict(cv_lasso, s = "lambda.1se", type = "nonzero") # get nonzero coefficients
if (nz_ind$lambda.1se[1]==1) {
  nz_ind <- nz_ind %>% slice(-1) # drop intercept
}
coef_names <- rownames(coef(cv_lasso,s="lambda.1se")[[1]])
nz_coef_names <- coef_names[nz_ind$lambda.1se]
lasso_parents <- check_coefs(vec = nz_coef_names, prefixes=colnames(dat_df[,-p_outcome]))
shds[2,2] <- bnl_shd(lasso_parents=lasso_parents, z_true = z_true)
### bnlearn
gs_fit <- learn.mb(dat_df, paste0("V",toString(p_outcome)), method="gs")
gs_parents <- as.numeric(colnames(dat_df[,-p_outcome])%in% gs_fit)
shds[2,3] <- bnl_shd(gs_parents, z_true)
iamb_fit <- learn.mb(dat_df, paste0("V",toString(p_outcome)), method="iamb")
iamb_parents <- as.numeric(colnames(dat_df[,-p_outcome])%in% iamb_fit)
shds[2,4] <- bnl_shd(iamb_parents, z_true)

### Fatigue
p_outcome <- which(colnames(dat)=="Fatigue")
z_true <- as.numeric( colnames(dat %>% select(-Fatigue)) %in% 
                   c("Lung_cancer","Coughing","Car_Accident", "Attention_Disorder")  )
## Gibbs
R <- 10000
B <- 200 # burn-in
epsilon <- c(0.32,0.32)
set.seed(1)
fit <- mala_edge_sampler_rcpp(R=R,
                              x_j=dat$Fatigue,
                              X_minus = dat %>% select(-Fatigue) %>% as.matrix(),
                              b_j=1, rho_j=2,
                              c_j=1, d_j=1,
                              step_size=epsilon,
                              t_init = rep(1,2),
                              stops=200)
colMeans(fit$t_accept[-(1:B),])
saveRDS(fit, "lucas_results/fatigue.rds")
### extract parent samples
parent_strings <- apply(fit$z[-(1:B),],1,paste,collapse="")
### MAP estimate
MAP <- as.numeric(strsplit(names(which.max(table(parent_strings))),split="")[[1]])
## log-posterior
log_posterior_fatigue <- log_posterior_samples(t=fit$t[-(1:B),],
                                                z=fit$z[-(1:B),],
                                                x_j = dat$Fatigue,
                                                X_minus = dat %>% select(-Fatigue) %>% as.matrix(),
                                                b_j=1, rho_j=2,
                                                c_j=1, d_j=1)
lp_df <- data.frame(iter = (B+1):R,
                    fatigue = log_posterior_fatigue)
trace_fatigue <- ggplot(lp_df, aes(x=iter,y=fatigue)) + geom_line() + 
  xlab("Iteration") + ylab("Log-Posterior") + labs(title = "Fatigue Log-Posterior MCMC Samples")  + theme_bw()
### posterior Hamming distance
#shd_samples <- apply(fit$z[-(1:B),],1,function(y) sum(y!=z_true))
shds[3,1] <- bnl_shd(MAP,z_true)
### lasso
formula_lasso <- as.formula(paste0("V", p_outcome, " ~ ."))
mm <- model.matrix(formula_lasso, data = dat_df)[,-1] 
cv_lasso <- cv.glmnet(x=mm, 
                      y=dat_df$V9, 
                      family = "multinomial",
                      type.multinomial="grouped",
                      alpha = 1) ## cross validation
coef(cv_lasso, s = "lambda.1se")
nz_ind <- predict(cv_lasso, s = "lambda.1se", type = "nonzero") # get nonzero coefficients
if (nz_ind$lambda.1se[1]==1) {
  nz_ind <- nz_ind %>% slice(-1) # drop intercept
}
coef_names <- rownames(coef(cv_lasso,s="lambda.1se")[[1]])
nz_coef_names <- coef_names[nz_ind$lambda.1se]
lasso_parents <- check_coefs(vec = nz_coef_names, prefixes=colnames(dat_df[,-p_outcome]))
shds[3,2] <- bnl_shd(lasso_parents=lasso_parents, z_true = z_true)
### bnlearn
gs_fit <- learn.mb(dat_df, paste0("V",toString(p_outcome)), method="gs")
gs_parents <- as.numeric(colnames(dat_df[,-p_outcome])%in% gs_fit)
shds[3,3] <- bnl_shd(gs_parents, z_true)
iamb_fit <- learn.mb(dat_df, paste0("V",toString(p_outcome)), method="iamb")
iamb_parents <- as.numeric(colnames(dat_df[,-p_outcome])%in% iamb_fit)
shds[3,4] <- bnl_shd(iamb_parents, z_true)

### Coughing
p_outcome <- which(colnames(dat)=="Coughing")
z_true <- as.numeric( colnames(dat %>% select(-Coughing)) %in% 
                        c("Allergy", "Lung_cancer", "Fatigue")  )
## Gibbs
R <- 10000
B <- 200 # burn-in
epsilon <- c(0.2,0.33)
set.seed(1)
fit <- mala_edge_sampler_rcpp(R=R,
                              x_j=dat$Coughing,
                              X_minus = dat %>% select(-Coughing) %>% as.matrix(),
                              b_j=1, rho_j=2,
                              c_j=1, d_j=1,
                              step_size=epsilon,
                              t_init = rep(1,2),
                              stops=200)
colMeans(fit$t_accept[-(1:B),])
saveRDS(fit, "lucas_results/coughing.rds")
### extract parent samples
parent_strings <- apply(fit$z[-(1:B),],1,paste,collapse="")
table(parent_strings)
### MAP estimate
MAP <- as.numeric(strsplit(names(which.max(table(parent_strings))),split="")[[1]])
## log-posterior
log_posterior_coughing <- log_posterior_samples(t=fit$t[-(1:B),],
                                               z=fit$z[-(1:B),],
                                               x_j = dat$Coughing,
                                               X_minus = dat %>% select(-Coughing) %>% as.matrix(),
                                               b_j=1, rho_j=2,
                                               c_j=1, d_j=1)
lp_df <- data.frame(iter = (B+1):R,
                    coughing = log_posterior_coughing)
trace_coughing <- ggplot(lp_df, aes(x=iter,y=coughing)) + geom_line() + 
  xlab("Iteration") + ylab("Log-Posterior") + labs(title = "Coughing Log-Posterior MCMC Samples")  + theme_bw()
### posterior Hamming distance
#shd_samples <- apply(fit$z[-(1:B),],1,function(y) sum(y!=z_true))
shds[4,1] <- bnl_shd(MAP,z_true)
### lasso
formula_lasso <- as.formula(paste0("V", p_outcome, " ~ ."))
mm <- model.matrix(formula_lasso, data = dat_df)[,-1] 
cv_lasso <- cv.glmnet(x=mm, 
                      y=dat_df$V11, 
                      family = "multinomial",
                      type.multinomial="grouped",
                      alpha = 1) ## cross validation
coef(cv_lasso, s = "lambda.1se")
nz_ind <- predict(cv_lasso, s = "lambda.1se", type = "nonzero") # get nonzero coefficients
if (nz_ind$lambda.1se[1]==1) {
  nz_ind <- nz_ind %>% slice(-1) # drop intercept
}
coef_names <- rownames(coef(cv_lasso,s="lambda.1se")[[1]])
nz_coef_names <- coef_names[nz_ind$lambda.1se]
lasso_parents <- check_coefs(vec = nz_coef_names, prefixes=colnames(dat_df[,-p_outcome]))
shds[4,2] <- bnl_shd(lasso_parents=lasso_parents, z_true = z_true)
### bnlearn
gs_fit <- learn.mb(dat_df, paste0("V",toString(p_outcome)), method="gs")
gs_parents <- as.numeric(colnames(dat_df[,-p_outcome])%in% gs_fit)
shds[4,3] <- bnl_shd(gs_parents, z_true)
iamb_fit <- learn.mb(dat_df, paste0("V",toString(p_outcome)), method="iamb")
iamb_parents <- as.numeric(colnames(dat_df[,-p_outcome])%in% iamb_fit)
shds[4,4] <- bnl_shd(iamb_parents, z_true)

### Allergy
p_outcome <- which(colnames(dat)=="Allergy")
z_true <- as.numeric( colnames(dat %>% select(-Allergy)) %in% 
                        c("Coughing", "Lung_cancer")  )
## Gibbs
R <- 10000
B <- 200 # burn-in
epsilon <- c(0.65,0.35)
set.seed(1)
fit <- mala_edge_sampler_rcpp(R=R,
                              x_j=dat$Allergy,
                              X_minus = dat %>% select(-Allergy) %>% as.matrix(),
                              b_j=1, rho_j=2,
                              c_j=1, d_j=1,
                              step_size=epsilon,
                              t_init = rep(1,2),
                              stops=200)
saveRDS(fit, "lucas_results/allergy.rds")
colMeans(fit$t_accept[-(1:B),])
### extract parent samples
parent_strings <- apply(fit$z[-(1:B),],1,paste,collapse="")
table(parent_strings)
### MAP estimate
MAP <- as.numeric(strsplit(names(which.max(table(parent_strings))),split="")[[1]])
## log-posterior
log_posterior_allergy <- log_posterior_samples(t=fit$t[-(1:B),],
                                                z=fit$z[-(1:B),],
                                                x_j = dat$Allergy,
                                                X_minus = dat %>% select(-Allergy) %>% as.matrix(),
                                                b_j=1, rho_j=2,
                                                c_j=1, d_j=1)
lp_df <- data.frame(iter = (B+1):R,
                    allergy = log_posterior_allergy)
trace_allergy <- ggplot(lp_df, aes(x=iter,y=allergy)) + geom_line() + 
  xlab("Iteration") + ylab("Log-Posterior") + labs(title = "Allergy Log-Posterior MCMC Samples")  + theme_bw()
### posterior Hamming distance
#shd_samples <- apply(fit$z[-(1:B),],1,function(y) sum(y!=z_true))
shds[5,1] <- bnl_shd(MAP,z_true)
### lasso
formula_lasso <- as.formula(paste0("V", p_outcome, " ~ ."))
mm <- model.matrix(formula_lasso, data = dat_df)[,-1] 
cv_lasso <- cv.glmnet(x=mm, 
                      y=dat_df$V10, 
                      family = "multinomial",
                      type.multinomial="grouped",
                      alpha = 1) ## cross validation
coef(cv_lasso, s = "lambda.1se")
nz_ind <- predict(cv_lasso, s = "lambda.1se", type = "nonzero") # get nonzero coefficients
if (nz_ind$lambda.1se[1]==1) {
  nz_ind <- nz_ind %>% slice(-1) # drop intercept
}
coef_names <- rownames(coef(cv_lasso,s="lambda.1se")[[1]])
nz_coef_names <- coef_names[nz_ind$lambda.1se]
lasso_parents <- check_coefs(vec = nz_coef_names, prefixes=colnames(dat_df[,-p_outcome]))
shds[5,2] <- bnl_shd(lasso_parents=lasso_parents, z_true = z_true)
### bnlearn
gs_fit <- learn.mb(dat_df, paste0("V",toString(p_outcome)), method="gs")
gs_parents <- as.numeric(colnames(dat_df[,-p_outcome])%in% gs_fit)
shds[5,3] <- bnl_shd(gs_parents, z_true)
iamb_fit <- learn.mb(dat_df, paste0("V",toString(p_outcome)), method="iamb")
iamb_parents <- as.numeric(colnames(dat_df[,-p_outcome])%in% iamb_fit)
shds[5,4] <- bnl_shd(iamb_parents, z_true)

### Lung Cancer
p_outcome <- which(colnames(dat)=="Lung_cancer")
z_true <- as.numeric( colnames(dat %>% select(-Lung_cancer)) %in% 
                        c("Smoking", "Genetics", "Fatigue", "Coughing", "Allergy")  )
## Gibbs
R <- 10000
B <- 200 # burn-in
epsilon <- c(0.11,0.15)
set.seed(1)
fit <- mala_edge_sampler_rcpp(R=R,
                              x_j=dat$Lung_cancer,
                              X_minus = dat %>% select(-Lung_cancer) %>% as.matrix(),
                              b_j=1, rho_j=2,
                              c_j=1, d_j=1,
                              step_size=epsilon,
                              t_init = rep(1,2),
                              stops=200)
colMeans(fit$t_accept[-(1:B),])
saveRDS(fit, "lucas_results/cancer.rds")
### extract parent samples
parent_strings <- apply(fit$z[-(1:B),],1,paste,collapse="")
table(parent_strings)
### MAP estimate
MAP <- as.numeric(strsplit(names(which.max(table(parent_strings))),split="")[[1]])
### posterior Hamming distance
#shd_samples <- apply(fit$z[-(1:B),],1,function(y) sum(y!=z_true))
shds[6,1] <- bnl_shd(MAP,z_true)
### edge probablities
edge_probs <- colMeans(fit$z[-(1:B),])
### traceplots
t_df <- data.frame(t1 = fit$t[-(1:B),1],
                   t2 = fit$t[-(1:B),2],
                   iter = (B+1):R)
t1_tp <- ggplot(t_df, aes(x=iter,y=t1)) + geom_line() + 
  xlab("Iteration") + ylab(TeX("$\\textbf{t}_j(1)$"))  + theme_bw()
t2_tp <- ggplot(t_df, aes(x=iter,y=t2)) + geom_line() +
  xlab("Iteration") + ylab(TeX("$\\textbf{t}_j(2)$")) + theme_bw()
plot_grid(t1_tp, t2_tp, nrow=1, labels = c("(a)", "(b)"))
coda::effectiveSize(fit$t[-(1:B),])
### log posterior
log_posterior_cancer <- log_posterior_samples(t=fit$t[-(1:B),],
                                               z=fit$z[-(1:B),],
                                               x_j = dat$Lung_cancer,
                                               X_minus = dat %>% select(-Lung_cancer) %>% as.matrix(),
                                               b_j=1, rho_j=2,
                                               c_j=1, d_j=1)
lp_df <- data.frame(iter = (B+1):R,
                    cancer = log_posterior_cancer)
trace_cancer <- ggplot(lp_df, aes(x=iter,y=cancer)) + geom_line() + 
  xlab("Iteration") + ylab("Log-Posterior") + labs(title = "Lung Cancer Log-Posterior MCMC Samples")  + theme_bw()
### lasso
formula_lasso <- as.formula(paste0("V", p_outcome, " ~ ."))
mm <- model.matrix(formula_lasso, data = dat_df)[,-1] 
cv_lasso <- cv.glmnet(x=mm, 
                      y=dat_df$V10, 
                      family = "multinomial",
                      type.multinomial="grouped",
                      alpha = 1) ## cross validation
coef(cv_lasso, s = "lambda.1se")
nz_ind <- predict(cv_lasso, s = "lambda.1se", type = "nonzero") # get nonzero coefficients
if (nz_ind$lambda.1se[1]==1) {
  nz_ind <- nz_ind %>% slice(-1) # drop intercept
}
coef_names <- rownames(coef(cv_lasso,s="lambda.1se")[[1]])
nz_coef_names <- coef_names[nz_ind$lambda.1se]
lasso_parents <- check_coefs(vec = nz_coef_names, prefixes=colnames(dat_df[,-p_outcome]))
shds[6,2] <- bnl_shd(lasso_parents=lasso_parents, z_true = z_true)
### bnlearn
gs_fit <- learn.mb(dat_df, paste0("V",toString(p_outcome)), method="gs")
gs_parents <- as.numeric(colnames(dat_df[,-p_outcome])%in% gs_fit)
shds[6,3] <- bnl_shd(gs_parents, z_true)
iamb_fit <- learn.mb(dat_df, paste0("V",toString(p_outcome)), method="iamb")
iamb_parents <- as.numeric(colnames(dat_df[,-p_outcome])%in% iamb_fit)
shds[6,4] <- bnl_shd(iamb_parents, z_true)
