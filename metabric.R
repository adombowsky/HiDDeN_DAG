# info: https://www.kaggle.com/datasets/raghadalharbi/breast-cancer-gene-expression-profiles-metabric
# create a folder "data" and save data as "METABRIC_RNA_Mutatation.csv"
source("r/mala_edges.R")
library(dplyr)
library(bnlearn)
library(reshape2)
cancer_data <- read.csv("data/METABRIC_RNA_Mutation.csv")
# remove sarcoma patient
cancer_data <- cancer_data[cancer_data$cancer_type=="Breast Cancer",]
# get variables
dat <- cancer_data %>%
  select(cancer_type_detailed, chemotherapy, type_of_breast_surgery, death_from_cancer, # outcomes
         age_at_diagnosis, tumor_size, radio_therapy, integrative_cluster,
         primary_tumor_laterality, inferred_menopausal_state, hormone_therapy,
         cohort)
# remove NAs
dat <- dat[-which(is.na(dat$tumor_size)),] # tumor size
dat <- dat[-which(dat$type_of_breast_surgery==""),] # type of breast surgery
dat <- dat[-which(dat$death_from_cancer==""),] # type of breast surgery
# dichotomize age
age_breaks <- c(20, 30, 40, 50, 60, 70, 80, 90, Inf)
dat$age_at_diagnosis <- cut(dat$age_at_diagnosis,breaks=age_breaks,right=F)
# dichotomize tumor size
tumor_size_breaks <- c(0, quantile(dat$tumor_size,probs=c(.25,.5,.75)), max(dat$tumor_size))
dat$tumor_size <- cut(dat$tumor_size,breaks=tumor_size_breaks,right=T)
# convert to HiDDeN format
dat <- dat %>% 
  mutate(cancer_type_detailed = as.integer(factor(cancer_type_detailed))-1,
         type_of_breast_surgery = as.integer(factor(type_of_breast_surgery))-1,
         death_from_cancer = as.integer(factor(death_from_cancer))-1,
         age_at_diagnosis = as.integer(age_at_diagnosis)-1,
         tumor_size = as.integer(tumor_size)-1,
         integrative_cluster = as.integer(factor(integrative_cluster))-1,
         primary_tumor_laterality = as.integer(factor(primary_tumor_laterality))-1,
         inferred_menopausal_state = as.integer(factor(inferred_menopausal_state))-1,
         cohort = cohort-1)

### preliminary function for lasso
extract_lasso_names <- function(variable_names) {
  cleaned_names <- gsub("\\d+$", "", variable_names)
  unique_names <- unique(cleaned_names)
  return(unique_names)
}


### fitting HiDDeN to the data ###

## Surgery
R <- 10000
B <- 500 # burn-in
k_j <- length(table(dat$type_of_breast_surgery))-1
epsilon <- c(0.175,0.25)
set.seed(1)
fit <- mala_edge_sampler_rcpp(R=R,
                              x_j=dat$type_of_breast_surgery,
                              X_minus = dat %>% select(cancer_type_detailed,
                                                       age_at_diagnosis,
                                                       tumor_size,
                                                       radio_therapy,
                                                       integrative_cluster,
                                                       primary_tumor_laterality,
                                                       inferred_menopausal_state,
                                                       hormone_therapy,
                                                       cohort) %>% as.matrix(),
                              b_j=1, rho_j=k_j+1,
                              c_j=1, d_j=1,
                              step_size=epsilon,
                              t_init = rep(1,k_j+1),
                              stops=200)
saveRDS(fit,"METABRIC_results/surgery.rds")
colMeans(fit$t_accept[-(1:B),])
parent_strings <- apply(fit$z[-(1:B),],1,paste,collapse="")
table(parent_strings)

### median probability model
edges_surgery <- colMeans(fit$z[-(1:B),])

### lasso
formula_lasso <- as.formula("type_of_breast_surgery ~ .")
mm <- model.matrix(formula_lasso, data = dat %>%
                     select(-death_from_cancer, -chemotherapy) %>%
                     mutate(across(.fns = as.factor)))[,-1] 
cv_lasso <- cv.glmnet(x=mm, 
                      y=as.factor(dat$type_of_breast_surgery), 
                      family = "multinomial",
                      type.multinomial="grouped",
                      alpha = 1) ## cross validation
coef_names <- rownames(coef(cv_lasso,s="lambda.1se")[[1]])
lasso_parents_surgery <- extract_lasso_names(coef_names)

### bnlearn
blacklist <- data.frame(
  from = rep("type_of_breast_surgery",9),
  to = colnames(dat %>% 
                  select(-type_of_breast_surgery, 
                         -chemotherapy, 
                         -death_from_cancer))
)
hc_fit <- hc(dat %>% select(-death_from_cancer, -chemotherapy) %>%
                    mutate(across(.fns = as.factor)),
             blacklist = blacklist)
tabu_fit <- tabu(dat %>% select(-death_from_cancer, -chemotherapy) %>%
               mutate(across(.fns = as.factor)),
             blacklist = blacklist)

## log-posterior 
log_posterior_surgery <- log_posterior_samples(t=fit$t[-(1:B),],
                                               z=fit$z[-(1:B),],
                                               x_j = dat$type_of_breast_surgery,
                                               X_minus = dat %>% select(cancer_type_detailed,
                                                                        age_at_diagnosis,
                                                                        tumor_size,
                                                                        radio_therapy,
                                                                        integrative_cluster,
                                                                        primary_tumor_laterality,
                                                                        inferred_menopausal_state,
                                                                        hormone_therapy,
                                                                        cohort) %>% as.matrix(),
                                               b_j=1, rho_j=k_j+1,
                                               c_j=1, d_j=1)

## Chemo
R <- 10000
B <- 500 # burn-in
k_j <- length(table(dat$chemotherapy))-1
epsilon <- c(.095,.045)
set.seed(1)
fit <- mala_edge_sampler_rcpp(R=R,
                              x_j=dat$chemotherapy,
                              X_minus = dat %>% select(cancer_type_detailed,
                                                       age_at_diagnosis,
                                                       tumor_size,
                                                       radio_therapy,
                                                       integrative_cluster,
                                                       primary_tumor_laterality,
                                                       inferred_menopausal_state,
                                                       hormone_therapy,
                                                       cohort,
                                                       type_of_breast_surgery) %>% as.matrix(),
                              b_j=1, rho_j=k_j+1,
                              c_j=1, d_j=1,
                              step_size=epsilon,
                              t_init = rep(1,k_j+1),
                              stops=200)
saveRDS(fit,"METABRIC_results/chemo.rds")
colMeans(fit$t_accept[-(1:B),])
parent_strings <- apply(fit$z[-(1:B),],1,paste,collapse="")
table(parent_strings)

### median probability model
edges_chemo <- colMeans(fit$z[-(1:B),])

### lasso
formula_lasso <- as.formula("chemotherapy ~ .")
mm <- model.matrix(formula_lasso, data = dat %>%
                     select(-death_from_cancer) %>%
                     mutate(across(.fns = as.factor)))[,-1] 
cv_lasso <- cv.glmnet(x=mm, 
                      y=as.factor(dat$chemotherapy), 
                      family = "multinomial",
                      type.multinomial="grouped",
                      alpha = 1) ## cross validation
coef_names <- rownames(coef(cv_lasso,s="lambda.1se")[[1]])
lasso_parents_chemo <- extract_lasso_names(coef_names)

### bnlearn
blacklist <- data.frame(
  from = rep("chemotherapy",10),
  to = colnames(dat %>% 
                  select(-chemotherapy, 
                         -death_from_cancer))
)
hc_fit <- hc(dat %>% select(-death_from_cancer) %>%
               mutate(across(.fns = as.factor)),
             blacklist = blacklist)
tabu_fit <- tabu(dat %>% select(-death_from_cancer) %>%
                   mutate(across(.fns = as.factor)),
                 blacklist = blacklist)

## log-posterior 
log_posterior_chemo <- log_posterior_samples(t=fit$t[-(1:B),],
                                               z=fit$z[-(1:B),],
                                               x_j = dat$chemotherapy,
                                               X_minus = dat %>% select(cancer_type_detailed,
                                                                        age_at_diagnosis,
                                                                        tumor_size,
                                                                        radio_therapy,
                                                                        integrative_cluster,
                                                                        primary_tumor_laterality,
                                                                        inferred_menopausal_state,
                                                                        hormone_therapy,
                                                                        cohort,
                                                                        type_of_breast_surgery) %>% as.matrix(),
                                               b_j=1, rho_j=k_j+1,
                                               c_j=1, d_j=1)

## Death
R <- 10000
B <- 500 # burn-in
k_j <- length(table(dat$death_from_cancer))-1
epsilon <- c(.35,.2,.35)
set.seed(1)
fit <- mala_edge_sampler_rcpp(R=R,
                              x_j=dat$death_from_cancer,
                              X_minus = dat %>% select(cancer_type_detailed,
                                                       age_at_diagnosis,
                                                       tumor_size,
                                                       radio_therapy,
                                                       integrative_cluster,
                                                       primary_tumor_laterality,
                                                       inferred_menopausal_state,
                                                       hormone_therapy,
                                                       cohort,
                                                       type_of_breast_surgery,
                                                       chemotherapy) %>% as.matrix(),
                              b_j=1, rho_j=k_j+1,
                              c_j=1, d_j=1,
                              step_size=epsilon,
                              t_init = rep(1,k_j+1),
                              stops=200)
saveRDS(fit,"METABRIC_results/death.rds")
colMeans(fit$t_accept[-(1:B),])
parent_strings <- apply(fit$z[-(1:B),],1,paste,collapse="")
table(parent_strings)

### median probability model
edges_death <- colMeans(fit$z[-(1:B),])

### lasso
formula_lasso <- as.formula("death_from_cancer ~ .")
mm <- model.matrix(formula_lasso, data = dat %>% mutate(across(.fns = as.factor)))[,-1] 
cv_lasso <- cv.glmnet(x=mm, 
                      y=as.factor(dat$death_from_cancer), 
                      family = "multinomial",
                      type.multinomial="grouped",
                      alpha = 1) ## cross validation
coef_names <- rownames(coef(cv_lasso,s="lambda.1se")[[1]])
lasso_parents_death <- extract_lasso_names(coef_names)

### bnlearn
blacklist <- data.frame(
  from = rep("death_from_cancer",11),
  to = colnames(dat %>% 
                  select(-death_from_cancer))
)
hc_fit <- hc(dat %>%
               mutate(across(.fns = as.factor)),
             blacklist = blacklist)
tabu_fit <- tabu(dat %>%
                   mutate(across(.fns = as.factor)),
                 blacklist = blacklist)

## log-posterior 
log_posterior_death <- log_posterior_samples(t=fit$t[-(1:B),],
                                             z=fit$z[-(1:B),],
                                             x_j = dat$death_from_cancer,
                                             X_minus = dat %>% select(cancer_type_detailed,
                                                                      age_at_diagnosis,
                                                                      tumor_size,
                                                                      radio_therapy,
                                                                      integrative_cluster,
                                                                      primary_tumor_laterality,
                                                                      inferred_menopausal_state,
                                                                      hormone_therapy,
                                                                      cohort,
                                                                      type_of_breast_surgery,
                                                                      chemotherapy) %>% as.matrix(),
                                             b_j=1, rho_j=k_j+1,
                                             c_j=1, d_j=1)

### plotting edge probabilities
edge_probs <- matrix(c(c(edges_surgery,0,0,0),
                c(edges_chemo,0,0),
                c(edges_death,0)),
                ncol=3)
rownames(edge_probs) = c("CTD", "AAD", "TSI", "RTH", "INC", "PTL", 
                         "IMS", "HTH", "COH", "TBS", "CHT", "DFC")
colnames(edge_probs) = c("TBS", "CHT", "DFC")
edge_probs <- edge_probs[nrow(edge_probs):1,] # reversing order for plotting
edges_long <- melt(edge_probs)
ggplot(edges_long,
  aes(x = Var2,y = Var1,fill = value)) +
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  geom_text(aes(label = sprintf("%.3f", value)), color = "black", size = 4) +
  scale_fill_gradient(
    low = "#F0F8FF", 
    high = "#F22105" 
  ) +
  labs(
    title = "Posterior edge probabilities",
    x = "Outcomes",
    y = "Variables",
    fill = "Prob"
  ) +
  # Ensure the tiles are square.
  coord_fixed() +
  # Customize the theme to make it look clean.
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

### plotting log posteriors
lp_df <- data.frame(iter = (B+1):R,
                    surgery = log_posterior_surgery,
                    chemo = log_posterior_chemo,
                    death = log_posterior_death)
trace_surgery <- ggplot(lp_df, aes(x=iter,y=surgery)) + geom_line() + 
  xlab("Iteration") + ylab("Log-Posterior") + labs(title = "Type of Breast Surgery")  + theme_bw()
trace_chemo <- ggplot(lp_df, aes(x=iter,y=chemo)) + geom_line() +
  xlab("Iteration") + ylab("Log-Posterior") + labs(title = "Chemotherapy")  + theme_bw()
trace_death <- ggplot(lp_df, aes(x=iter,y=death)) + geom_line() +
  xlab("Iteration") + ylab("Log-Posterior") + labs(title = "Death from Cancer")  + theme_bw()
plot_grid(trace_surgery, trace_chemo, trace_death, ncol=1, labels = c("(a)", "(b)", "(c)"))
