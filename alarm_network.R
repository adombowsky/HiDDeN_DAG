# preliminaries
library(bnlearn)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer) # for colors
source("r/mala_graph.R")

# MCMC samples
R <- 50000
B <- 200

# loading in data
full_data <- bnlearn::alarm
dat <- full_data %>%
  dplyr::select(LVF, HYP, HIST, LVV, STKV, PCWP, CVP)
K <- apply(dat,2,function(x) length(table(x)))-1

# subsetting the data
seed <- 1 # random number seed
n <- 200 # sample size, < 20000
set.seed(seed)
ids <- sample(1:nrow(dat), size = n, replace = F)
x <- dat[ids,]

# calibrate to be numeric
LNH_to_numeric <- function(y) {
  z <- rep(0,length(y))
  for (i in 1:length(y)) {
    if (y[i]=="LOW") {
      z[i] <- 1
    } else if (y[i]=="NORMAL") {
      z[i] <- 2
    } else {
      z[i] <- 3
    }
  }
  return(z)
}
x <- x %>%
  mutate(LVF = as.numeric(LVF),
         HYP = as.numeric(HYP),
         HIST = as.numeric(HIST),
         LVV = LNH_to_numeric(LVV),
         STKV = LNH_to_numeric(STKV),
         PCWP = LNH_to_numeric(PCWP),
         CVP = LNH_to_numeric(CVP))

### focus on intermediary variables, test for edges to parents
# HIST
## test 3 graphs
# under graph 1
mct_hist1 <- as.matrix(as.vector(table(x$LVF)))
ct_hist1 <- as.matrix(table(x$HIST,x$LVF))
par_grid_hist1 <- expand.grid(list(0:K["LVF"]))
par_child_grid_hist1 <- expand.grid(list(0:K["HIST"],0:K["LVF"]))

# under graph 2
mct_hist2 <- as.matrix(as.vector(table(x$HYP)))
ct_hist2 <- as.matrix(table(x$HIST,x$HYP))
par_grid_hist2 <- expand.grid(list(0:K["HYP"]))
par_child_grid_hist2 <- expand.grid(list(0:K["HIST"],0:K["HYP"]))

# under graph 3
mct_hist3 <- as.matrix(as.vector(table(x$LVF, x$HYP)))
ct_hist3 <- as.matrix(as.vector(table(x$HIST,x$LVF, x$HYP)))
par_grid_hist3 <- expand.grid(list(0:K["LVF"], 0:K["HYP"]))
par_child_grid_hist3 <- expand.grid(list(0:K["HIST"],0:K["LVF"],0:K["HYP"]))

step_size <- 3*10^(-1)*rep(1,K["HIST"]+1)
fit_hist <- mala_graph_sampler(R=R,
                          ct=list(ct_hist1,ct_hist2,ct_hist3),
                          mct = list(mct_hist1, mct_hist2, mct_hist3),
                          par_grid = list(par_grid_hist1, par_grid_hist2, par_grid_hist3),
                          par_child_grid = list(par_child_grid_hist1, par_child_grid_hist2, par_child_grid_hist3),
                          b=1,rho=K["HIST"]+1,k=K["HIST"],
                          prior_probs = c(1/3,1/3,1/3),
                          step_size = step_size)
pset_hist <- fit_hist$pset[-(1:B)]

# LVV
## test 7 graphs
# under graph 1
mct_lvv1 <- as.matrix(as.vector(table(x$LVF)))
ct_lvv1 <- as.matrix(table(x$LVV,x$LVF))
par_grid_lvv1 <- expand.grid(list(0:K["LVF"]))
par_child_grid_lvv1 <- expand.grid(list(0:K["LVV"],0:K["LVF"]))

# under graph 2
mct_lvv2 <- as.matrix(as.vector(table(x$HYP)))
ct_lvv2 <- as.matrix(table(x$LVV,x$HYP))
par_grid_lvv2 <- expand.grid(list(0:K["HYP"]))
par_child_grid_lvv2 <- expand.grid(list(0:K["LVV"],0:K["HYP"]))

# under graph 3
mct_lvv3 <- as.matrix(as.vector(table(x$HIST)))
ct_lvv3 <- as.matrix(table(x$LVV,x$HIST))
par_grid_lvv3 <- expand.grid(list(0:K["HIST"]))
par_child_grid_lvv3 <- expand.grid(list(0:K["LVV"],0:K["HIST"]))

# under graph 4
mct_lvv4 <- as.matrix(as.vector(table(x$LVF, x$HYP)))
ct_lvv4 <- as.matrix(as.vector(table(x$LVV,x$LVF, x$HYP)))
par_grid_lvv4 <- expand.grid(list(0:K["LVF"], 0:K["HYP"]))
par_child_grid_lvv4 <- expand.grid(list(0:K["LVV"],0:K["LVF"],0:K["HYP"]))

# under graph 5
mct_lvv5 <- as.matrix(as.vector(table(x$LVF, x$HIST)))
ct_lvv5 <- as.matrix(as.vector(table(x$LVV,x$LVF, x$HIST)))
par_grid_lvv5 <- expand.grid(list(0:K["LVF"], 0:K["HIST"]))
par_child_grid_lvv5 <- expand.grid(list(0:K["LVV"],0:K["LVF"],0:K["HIST"]))

# under graph 6
mct_lvv6 <- as.matrix(as.vector(table(x$HIST, x$HYP)))
ct_lvv6 <- as.matrix(as.vector(table(x$LVV,x$HIST, x$HYP)))
par_grid_lvv6 <- expand.grid(list(0:K["HIST"], 0:K["HYP"]))
par_child_grid_lvv6 <- expand.grid(list(0:K["LVV"],0:K["HIST"],0:K["HYP"]))

# under graph 7
mct_lvv7 <- as.matrix(as.vector(table(x$HIST, x$LVF, x$HYP)))
ct_lvv7 <- as.matrix(as.vector(table(x$LVV,x$HIST, x$LVF, x$HYP)))
par_grid_lvv7 <- expand.grid(list(0:K["HIST"], 0:K["LVF"], 0:K["HYP"]))
par_child_grid_lvv7 <- expand.grid(list(0:K["LVV"],0:K["HIST"], 0:K["LVF"], 0:K["HYP"]))

step_size <- 2*10^(-1)*rep(1,K["LVV"]+1)
fit_lvv <- mala_graph_sampler(R=R,
                               ct=list(ct_lvv1,ct_lvv2,ct_lvv3,ct_lvv4,ct_lvv5,ct_lvv6,ct_lvv7),
                               mct = list(mct_lvv1, mct_lvv2, mct_lvv3,mct_lvv4, mct_lvv5, mct_lvv6, mct_lvv7),
                               par_grid = list(par_grid_lvv1, par_grid_lvv2, par_grid_lvv3,
                                               par_grid_lvv4, par_grid_lvv5, par_grid_lvv6,
                                               par_grid_lvv7),
                               par_child_grid = list(par_child_grid_lvv1, par_child_grid_lvv2, par_child_grid_lvv3,
                                                     par_child_grid_lvv4, par_child_grid_lvv5, par_child_grid_lvv6,
                                                     par_child_grid_lvv7),
                               b=1,rho=K["LVV"]+1,k=K["LVV"],
                               prior_probs = rep(1/7,7),
                               step_size = step_size)
colMeans(fit_lvv$t_accept)
pset_lvv <- fit_lvv$pset[-(1:B)]

# STKV
## test 3 graphs
# under graph 1
mct_stkv1 <- as.matrix(as.vector(table(x$LVF)))
ct_stkv1 <- as.matrix(table(x$STKV,x$LVF))
par_grid_stkv1 <- expand.grid(list(0:K["LVF"]))
par_child_grid_stkv1 <- expand.grid(list(0:K["STKV"],0:K["LVF"]))

# under graph 2
mct_stkv2 <- as.matrix(as.vector(table(x$HYP)))
ct_stkv2 <- as.matrix(table(x$STKV,x$HYP))
par_grid_stkv2 <- expand.grid(list(0:K["HYP"]))
par_child_grid_stkv2 <- expand.grid(list(0:K["STKV"],0:K["HYP"]))

# under graph 3
mct_stkv3 <- as.matrix(as.vector(table(x$HIST)))
ct_stkv3 <- as.matrix(table(x$STKV,x$HIST))
par_grid_stkv3 <- expand.grid(list(0:K["HIST"]))
par_child_grid_stkv3 <- expand.grid(list(0:K["STKV"],0:K["HIST"]))

# under graph 4
mct_stkv4 <- as.matrix(as.vector(table(x$LVV)))
ct_stkv4 <- as.matrix(table(x$STKV,x$LVV))
par_grid_stkv4 <- expand.grid(list(0:K["LVV"]))
par_child_grid_stkv4 <- expand.grid(list(0:K["STKV"],0:K["LVV"]))

# under graph 5
mct_stkv5 <- as.matrix(as.vector(table(x$LVF, x$HYP)))
ct_stkv5 <- as.matrix(as.vector(table(x$STKV,x$LVF, x$HYP)))
par_grid_stkv5 <- expand.grid(list(0:K["LVF"], 0:K["HYP"]))
par_child_grid_stkv5 <- expand.grid(list(0:K["STKV"],0:K["LVF"],0:K["HYP"]))

# under graph 6
mct_stkv6 <- as.matrix(as.vector(table(x$LVF, x$HIST)))
ct_stkv6 <- as.matrix(as.vector(table(x$STKV,x$LVF, x$HIST)))
par_grid_stkv6 <- expand.grid(list(0:K["LVF"], 0:K["HIST"]))
par_child_grid_stkv6 <- expand.grid(list(0:K["STKV"],0:K["LVF"],0:K["HIST"]))

# under graph 7
mct_stkv7 <- as.matrix(as.vector(table(x$HYP, x$HIST)))
ct_stkv7 <- as.matrix(as.vector(table(x$STKV,x$HYP, x$HIST)))
par_grid_stkv7 <- expand.grid(list(0:K["HYP"], 0:K["HIST"]))
par_child_grid_stkv7 <- expand.grid(list(0:K["STKV"],0:K["HYP"],0:K["HIST"]))

# under graph 8
mct_stkv8 <- as.matrix(as.vector(table(x$HYP, x$LVV)))
ct_stkv8 <- as.matrix(as.vector(table(x$STKV,x$HYP, x$LVV)))
par_grid_stkv8 <- expand.grid(list(0:K["HYP"], 0:K["LVV"]))
par_child_grid_stkv8 <- expand.grid(list(0:K["STKV"],0:K["HYP"],0:K["LVV"]))

# under graph 9
mct_stkv9 <- as.matrix(as.vector(table(x$LVV, x$HIST)))
ct_stkv9 <- as.matrix(as.vector(table(x$STKV,x$LVV, x$HIST)))
par_grid_stkv9 <- expand.grid(list(0:K["LVV"], 0:K["HIST"]))
par_child_grid_stkv9 <- expand.grid(list(0:K["STKV"],0:K["LVV"],0:K["HIST"]))

# under graph 10
mct_stkv10 <- as.matrix(as.vector(table(x$LVF, x$LVV)))
ct_stkv10 <- as.matrix(as.vector(table(x$STKV,x$LVF, x$LVV)))
par_grid_stkv10 <- expand.grid(list(0:K["LVF"], 0:K["LVV"]))
par_child_grid_stkv10 <- expand.grid(list(0:K["STKV"],0:K["LVF"],0:K["LVV"]))

# under graph 11
mct_stkv11 <- as.matrix(as.vector(table(x$HIST, x$LVF, x$HYP)))
ct_stkv11 <- as.matrix(as.vector(table(x$STKV,x$HIST, x$LVF, x$HYP)))
par_grid_stkv11 <- expand.grid(list(0:K["HIST"], 0:K["LVF"], 0:K["HYP"]))
par_child_grid_stkv11 <- expand.grid(list(0:K["STKV"],0:K["HIST"], 0:K["LVF"], 0:K["HYP"]))

# under graph 12
mct_stkv12 <- as.matrix(as.vector(table(x$LVV, x$LVF, x$HYP)))
ct_stkv12 <- as.matrix(as.vector(table(x$STKV,x$LVV, x$LVF, x$HYP)))
par_grid_stkv12 <- expand.grid(list(0:K["LVV"], 0:K["LVF"], 0:K["HYP"]))
par_child_grid_stkv12 <- expand.grid(list(0:K["STKV"],0:K["LVV"], 0:K["LVF"], 0:K["HYP"]))

# under graph 13
mct_stkv13 <- as.matrix(as.vector(table(x$HIST, x$LVV, x$HYP)))
ct_stkv13 <- as.matrix(as.vector(table(x$STKV,x$HIST, x$LVV, x$HYP)))
par_grid_stkv13 <- expand.grid(list(0:K["HIST"], 0:K["LVV"], 0:K["HYP"]))
par_child_grid_stkv13 <- expand.grid(list(0:K["STKV"],0:K["HIST"], 0:K["LVV"], 0:K["HYP"]))

# under graph 14
mct_stkv14 <- as.matrix(as.vector(table(x$HIST, x$LVF, x$LVV)))
ct_stkv14 <- as.matrix(as.vector(table(x$STKV,x$HIST, x$LVF, x$LVV)))
par_grid_stkv14 <- expand.grid(list(0:K["HIST"], 0:K["LVF"], 0:K["LVV"]))
par_child_grid_stkv14 <- expand.grid(list(0:K["STKV"],0:K["HIST"], 0:K["LVF"], 0:K["LVV"]))

# under graph 15
mct_stkv15 <- as.matrix(as.vector(table(x$HIST, x$LVF, x$HYP, x$LVV)))
ct_stkv15 <- as.matrix(as.vector(table(x$STKV,x$HIST, x$LVF, x$HYP, x$LVV)))
par_grid_stkv15 <- expand.grid(list(0:K["HIST"], 0:K["LVF"], 0:K["HYP"], 0:K["LVV"]))
par_child_grid_stkv15 <- expand.grid(list(0:K["STKV"],0:K["HIST"], 0:K["LVF"], 0:K["HYP"], 0:K["LVV"]))


step_size <- c(4*10^(-1), 4*10^(-1) ,2*10^(-1))
fit_stkv <- mala_graph_sampler(R=R,
                              ct=list(ct_stkv1,ct_stkv2,ct_stkv3,
                                      ct_stkv4,ct_stkv5,ct_stkv6,
                                      ct_stkv7,ct_stkv8,ct_stkv9,
                                      ct_stkv10,ct_stkv11,ct_stkv12,
                                      ct_stkv13,ct_stkv14,ct_stkv15),
                              mct = list(mct_stkv1, mct_stkv2, mct_stkv3,
                                         mct_stkv4, mct_stkv5, mct_stkv6,
                                         mct_stkv7, mct_stkv8, mct_stkv9,
                                         mct_stkv10, mct_stkv11, mct_stkv12,
                                         mct_stkv13, mct_stkv14, mct_stkv15),
                              par_grid = list(par_grid_stkv1, par_grid_stkv2, par_grid_stkv3,
                                              par_grid_stkv4, par_grid_stkv5, par_grid_stkv6,
                                              par_grid_stkv7, par_grid_stkv8, par_grid_stkv9,
                                              par_grid_stkv10, par_grid_stkv11, par_grid_stkv12,
                                              par_grid_stkv13, par_grid_stkv14, par_grid_stkv15),
                              par_child_grid = list(par_child_grid_stkv1, par_child_grid_stkv2, par_child_grid_stkv3,
                                                    par_child_grid_stkv4, par_child_grid_stkv5, par_child_grid_stkv6,
                                                    par_child_grid_stkv7, par_child_grid_stkv8, par_child_grid_stkv9,
                                                    par_child_grid_stkv10, par_child_grid_stkv11, par_child_grid_stkv12,
                                                    par_child_grid_stkv13, par_child_grid_stkv14, par_child_grid_stkv15),
                              b=1,rho=K["STKV"]+1,k=K["STKV"],
                              prior_probs = rep(1/15,15),
                              step_size = step_size)
colMeans(fit_stkv$t_accept) # acceptance probability
pset_stkv <- fit_stkv$pset[-(1:B)]


# comparing graphs
phat_hist <- round(table(factor(pset_hist, levels = 1:3))/length(pset_hist),3)
phat_lvv <- round(table(factor(pset_lvv,levels=1:7))/length(pset_lvv),3)
phat_stkv <- round(table(factor(pset_stkv,levels=1:15))/length(pset_stkv),3)

# for stkv, computing edge probabilities
tab_stkv_full <- table(factor(pset_stkv, levels = 1:15))
round(sum(tab_stkv_full[c(1,5,6,10,11,12,14,15)])/length(pset_stkv),3)
round(sum(tab_stkv_full[c(2, 5, 7, 8, 11, 12, 13,15)])/length(pset_stkv),3)

# computing edge probabilities

# comparison to score based techniques using bnlearn
# g_1 <- "[x_1][x_2][x_3|x_1]"
# g_1_dag <- model2network(g_1)
# g_2 <- "[x_1][x_2][x_3|x_2]"
# g_2_dag <- model2network(g_2)
# g_3 <- "[x_1][x_2][x_3|x_1:x_2]"
# g_3_dag <- model2network(g_3)
# ## BIC
# bnlearn::score(x=g_1_dag, data=data.frame(x_1 = as.factor(x$LVF),
#                                           x_2 = as.factor(x$HYP),
#                                           x_3 = as.factor(x$HIST)), type="bde")
# bnlearn::score(x=g_2_dag, data=data.frame(x_1 = as.factor(x$LVF),
#                                           x_2 = as.factor(x$HYP),
#                                           x_3 = as.factor(x$HIST)), type="bde")
# bnlearn::score(x=g_3_dag, data=data.frame(x_1 = as.factor(x$LVF),
#                                           x_2 = as.factor(x$HYP),
#                                           x_3 = as.factor(x$HIST)), type="bde")

