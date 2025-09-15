# preliminaries
library(bnlearn)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer) # for colors
source("r/mala_within_gibbs.R")

# MCMC samples
R <- 50000

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

# HIST
par_grid_hist <- expand.grid(list(0:K["LVF"]))
par_child_grid_hist <- expand.grid(list(0:K["HIST"],0:K["LVF"]))
mct_hist <- as.matrix(table(x$LVF))
ct_hist <- as.matrix(as.vector(table(x$HIST,x$LVF)))
step_size <- 3*10^(-1)*rep(1,K["HIST"]+1)
mala_hist <- mala_within_gibbs_sampler(R=R,
                                       ct=ct_hist,
                                       mct=mct_hist,
                                       par_grid=par_grid_hist,
                                       par_child_grid = par_child_grid_hist,
                                       b=1,rho=K["HIST"]+1,k=K["HIST"],
                                       step_size=step_size)
colMeans(mala_hist$t_accept)
samps_hist <- classifier_MCMC_samples(t_samps=mala_hist$t, 
                                      ct=ct_hist, 
                                      mct=mct_hist, 
                                      par_child_grid = par_child_grid_hist, 
                                      par_grid = par_grid_hist)
posterior_mean_hist <- colMeans(samps_hist)


# LVV
par_grid_lvv <- expand.grid(list(0:K["LVF"], 0:K["HYP"]))
par_child_grid_lvv <- expand.grid(list(0:K["LVV"],0:K["LVF"], 0:K["HYP"]))
mct_lvv <- as.matrix(table(x$LVF, x$HYP))
ct_lvv <- as.matrix(as.vector(table(x$LVV,x$LVF, x$HYP)))
step_size <- 2*10^(-1)*rep(1,K["LVV"]+1)
mala_lvv <- mala_within_gibbs_sampler(R=R,
                                       ct=ct_lvv,
                                       mct=mct_lvv,
                                       par_grid=par_grid_lvv,
                                       par_child_grid = par_child_grid_lvv,
                                       b=1,rho=K["LVV"]+1,k=K["LVV"],
                                       step_size=step_size)
colMeans(mala_lvv$t_accept)
samps_lvv <- classifier_MCMC_samples(t_samps=mala_lvv$t, 
                                      ct=ct_lvv, 
                                      mct=mct_lvv, 
                                      par_child_grid = par_child_grid_lvv, 
                                      par_grid = par_grid_lvv)
posterior_mean_lvv <- colMeans(samps_lvv)

# STKV
par_grid_stkv <- expand.grid(list(0:K["LVF"], 0:K["HYP"]))
par_child_grid_stkv <- expand.grid(list(0:K["STKV"],0:K["LVF"], 0:K["HYP"]))
mct_stkv <- as.matrix(table(x$LVF, x$HYP))
ct_stkv <- as.matrix(as.vector(table(x$STKV,x$LVF, x$HYP)))
step_size <- c(4*10^(-1), 4*10^(-1) ,2*10^(-1))
mala_stkv <- mala_within_gibbs_sampler(R=R,
                                      ct=ct_stkv,
                                      mct=mct_stkv,
                                      par_grid=par_grid_stkv,
                                      par_child_grid = par_child_grid_stkv,
                                      b=1,rho=K["STKV"]+1,k=K["STKV"],
                                      step_size=step_size)
colMeans(mala_stkv$t_accept)
samps_stkv <- classifier_MCMC_samples(t_samps=mala_stkv$t, 
                                     ct=ct_stkv, 
                                     mct=mct_stkv, 
                                     par_child_grid = par_child_grid_stkv, 
                                     par_grid = par_grid_stkv)
posterior_mean_stkv <- colMeans(samps_stkv)

# PCWP
par_grid_pcwp <- expand.grid(list(0:K["LVV"]))
par_child_grid_pcwp <- expand.grid(list(0:K["PCWP"],0:K["LVV"]))
mct_pcwp <- as.matrix(table(x$LVV))
ct_pcwp <- as.matrix(as.vector(table(x$PCWP,x$LVV)))
step_size <- 2*10^(-1)*rep(1,K["PCWP"]+1)
mala_pcwp <- mala_within_gibbs_sampler(R=R,
                                       ct=ct_pcwp,
                                       mct=mct_pcwp,
                                       par_grid=par_grid_pcwp,
                                       par_child_grid = par_child_grid_pcwp,
                                       b=1,rho=K["PCWP"]+1,k=K["PCWP"],
                                       step_size=step_size)
colMeans(mala_pcwp$t_accept)
samps_pcwp <- classifier_MCMC_samples(t_samps=mala_pcwp$t, 
                                      ct=ct_pcwp, 
                                      mct=mct_pcwp, 
                                      par_child_grid = par_child_grid_pcwp, 
                                      par_grid = par_grid_pcwp)
posterior_mean_pcwp <- colMeans(samps_pcwp)

# CVP
par_grid_cvp <- expand.grid(list(0:K["LVV"]))
par_child_grid_cvp <- expand.grid(list(0:K["CVP"],0:K["LVV"]))
mct_cvp <- as.matrix(table(x$LVV))
ct_cvp <- as.matrix(as.vector(table(x$CVP,x$LVV)))
step_size <- c(0.2,0.25,0.15)
mala_cvp <- mala_within_gibbs_sampler(R=R,
                                       ct=ct_cvp,
                                       mct=mct_cvp,
                                       par_grid=par_grid_cvp,
                                       par_child_grid = par_child_grid_cvp,
                                       b=1,rho=K["CVP"]+1,k=K["CVP"],
                                       step_size=step_size)
colMeans(mala_cvp$t_accept)
samps_cvp <- classifier_MCMC_samples(t_samps=mala_cvp$t, 
                                      ct=ct_cvp, 
                                      mct=mct_cvp, 
                                      par_child_grid = par_child_grid_cvp, 
                                      par_grid = par_grid_cvp)
posterior_mean_cvp <- colMeans(samps_cvp)

# saving output
samps <- list(HIST = samps_hist,
              LVV = samps_lvv,
              STKV = samps_stkv,
              PCWP = samps_pcwp,
              CVP = samps_cvp)
posterior_means <- list(HIST=posterior_mean_hist,
                        LVV = posterior_mean_lvv,
                        STKV = posterior_mean_stkv,
                        PCWP = posterior_mean_pcwp,
                        CVP=posterior_mean_cvp)


### Plots ###

# HIST
plot_data <- data.frame()
posterior_mean_hist_mat <- matrix(posterior_mean_hist,ncol=K["HIST"]+1,byrow=T)
for (i in 1:nrow(par_grid_hist)) {
  current_pa <- par_grid_hist[i,]
  
  # Create a temporary data frame for this parent configuration
  temp_df <- data.frame(
    ChildNode = as.factor(0:(K["HIST"])),
    Pa = current_pa,
    Probability = sqrt(posterior_mean_hist_mat[i,])
  )
  
  # Append to the main plot data frame
  plot_data <- rbind(plot_data, temp_df)
}

plot_data <- plot_data %>%
  mutate(ParentConfig = paste0("LVF=", Pa))

bar_hist <- ggplot(plot_data, aes(x = ChildNode, y = Probability, fill = ParentConfig)) +
  geom_bar(stat = "identity", position = "dodge", color = "white", linewidth = 0.5) + # 'dodge' positions bars side-by-side
  labs(
    title = "HIST",
    x = " ",
    y = "Sqrt Cond Prob",
    fill = "Parent Values" # Legend title
  ) + scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + # Ensure y-axis from 0 to 1
  theme_bw() + # A clean theme
  theme(
    plot.title = element_text(hjust = 0.5,size=20), # Center and bold title
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(), # Remove vertical grid lines
    panel.grid.minor.x = element_blank()
  )

# LVV
plot_data <- data.frame()
posterior_mean_lvv_mat <- matrix(posterior_mean_lvv,ncol=K["LVV"]+1,byrow=T)
for (i in 1:nrow(par_grid_lvv)) {
  #current_pa <- par_grid_lvv[i,]
  
  # Create a temporary data frame for this parent configuration
  temp_df <- data.frame(
    ChildNode = 1:(K["LVV"]+1),
    Pa1 = par_grid_lvv[i,1],
    Pa2 = par_grid_lvv[i,2],
    Probability = sqrt(posterior_mean_lvv_mat[i,])
  )
  
  # Append to the main plot data frame
  plot_data <- rbind(plot_data, temp_df)
}

plot_data <- plot_data %>%
  mutate(ParentConfig = paste0("LVF=", Pa1, ", HYP=", Pa2))

bar_lvv <- ggplot(plot_data, aes(x = ChildNode, y = Probability, fill = ParentConfig)) +
  geom_bar(stat = "identity", position = "dodge", color = "white", linewidth = 0.5) + # 'dodge' positions bars side-by-side
  labs(
    title = "LVV",
    x = " ",
    y = "Sqrt Cond Prob",
    fill = "Parent Values" # Legend title
  ) + scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + # Ensure y-axis from 0 to 1
  theme_bw() + # A clean theme
  theme(
    plot.title = element_text(hjust = 0.5,size=20), # Center and bold title
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(), # Remove vertical grid lines
    panel.grid.minor.x = element_blank()
  )

# STKV
plot_data <- data.frame()
posterior_mean_stkv_mat <- matrix(posterior_mean_stkv,ncol=K["STKV"]+1,byrow=T)
for (i in 1:nrow(par_grid_stkv)) {
  #current_pa <- par_grid_lvv[i,]
  
  # Create a temporary data frame for this parent configuration
  temp_df <- data.frame(
    ChildNode = 1:(K["STKV"]+1),
    Pa1 = par_grid_stkv[i,1],
    Pa2 = par_grid_stkv[i,2],
    Probability = sqrt(posterior_mean_stkv_mat[i,])
  )
  
  # Append to the main plot data frame
  plot_data <- rbind(plot_data, temp_df)
}

plot_data <- plot_data %>%
  mutate(ParentConfig = paste0("LVF=", Pa1, ", HYP=", Pa2))

bar_stkv <- ggplot(plot_data, aes(x = ChildNode, y = Probability, fill = ParentConfig)) +
  geom_bar(stat = "identity", position = "dodge", color = "white", linewidth = 0.5) + # 'dodge' positions bars side-by-side
  labs(
    title = "STKV",
    x = " ",
    y = "Sqrt Cond Prob",
    fill = "Parent Values" # Legend title
  ) +  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + # Ensure y-axis from 0 to 1
  theme_bw() + # A clean theme
  theme(
    plot.title = element_text(hjust = 0.5,size=20), # Center and bold title
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(), # Remove vertical grid lines
    panel.grid.minor.x = element_blank()
  )

# PCWP
plot_data <- data.frame()
posterior_mean_pcwp_mat <- matrix(posterior_mean_pcwp,ncol=K["PCWP"]+1,byrow=T)
for (i in 1:nrow(par_grid_pcwp)) {
  current_pa <- par_grid_pcwp[i,]
  
  # Create a temporary data frame for this parent configuration
  temp_df <- data.frame(
    ChildNode = 1:(K["PCWP"]+1),
    Pa = current_pa,
    Probability = sqrt(posterior_mean_pcwp_mat[i,])
  )
  
  # Append to the main plot data frame
  plot_data <- rbind(plot_data, temp_df)
}

plot_data <- plot_data %>%
  mutate(ParentConfig = paste0("LVV=", Pa+1))

bar_pcwp <- ggplot(plot_data, aes(x = ChildNode, y = Probability, fill = ParentConfig)) +
  geom_bar(stat = "identity", position = "dodge", color = "white", linewidth = 0.5) + # 'dodge' positions bars side-by-side
  labs(
    title = "PCWP",
    x = " ",
    y = "Sqrt Cond Prob",
    fill = "Parent Values" # Legend title
  ) +  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + # Ensure y-axis from 0 to 1
  theme_bw() + # A clean theme
  theme(
    plot.title = element_text(hjust = 0.5,size=20), # Center and bold title
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(), # Remove vertical grid lines
    panel.grid.minor.x = element_blank()
  )


# CWP
plot_data <- data.frame()
posterior_mean_cvp_mat <- matrix(posterior_mean_cvp,ncol=K["CVP"]+1,byrow=T)
for (i in 1:nrow(par_grid_cvp)) {
  current_pa <- par_grid_cvp[i,]
  
  # Create a temporary data frame for this parent configuration
  temp_df <- data.frame(
    ChildNode = 1:(K["CVP"]+1),
    Pa = current_pa,
    Probability = sqrt(posterior_mean_cvp_mat[i,])
  )
  
  # Append to the main plot data frame
  plot_data <- rbind(plot_data, temp_df)
}

plot_data <- plot_data %>%
  mutate(ParentConfig = paste0("LVV=", Pa+1))

bar_cvp <- ggplot(plot_data, aes(x = ChildNode, y = Probability, fill = ParentConfig)) +
  geom_bar(stat = "identity", position = "dodge", color = "white", linewidth = 0.5) + # 'dodge' positions bars side-by-side
  labs(
    title = "CVP",
    x = " ",
    y = "Sqrt Cond Prob",
    fill = "Parent Values" # Legend title
  ) +  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + # Ensure y-axis from 0 to 1
  theme_bw() + # A clean theme
  theme(
    plot.title = element_text(hjust = 0.5,size=20), # Center and bold title
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(), # Remove vertical grid lines
    panel.grid.minor.x = element_blank()
  )

# plotting
plot_grid(bar_hist,
          bar_lvv,
          bar_stkv,
          bar_pcwp,
          bar_cvp,
          ncol=1)
