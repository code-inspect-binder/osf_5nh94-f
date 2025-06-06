library(ggplot2)
library(dplyr)
library(cowplot)
library(multibergm)
theme_set(theme_cowplot(font_size=12))
set.seed(42)

thresh <- "constK3" # vary between constK3, constK5, propK3, propK5
constraints <- if (substr(thresh,1,1) == "p")  ~edges else ~.

stat_labels <- c("edges","gwesp","gwnsp")
if (substr(thresh,1,1) == 'p'){
  stat_labels <- stat_labels[2:3]
}

group_fit <- readRDS(file.path("Output", thresh, "combined.RDS"))

# Remove burnIn iterations and apply thinning (default: no thinning)
burn_in        <- 1000
thin          <- 1
post_iters     <- seq(burn_in + 1, group_fit$mainIters, thin)
group_fit$params <- lapply(group_fit$params,
                          function(x) abind::asub(x, post_iters, 1))
n_stats <- length(group_fit$control$model$coef.names)
n_iters <- dim(group_fit$params$muPop)[1]

n_nets    <- length(group_fit$networks)
n_sample  <- 1000
aux_iters <- 100000

graphs <- lapply(group_fit$networks,
                 function(x) igraph::graph_from_adjacency_matrix(as.matrix(x)))
local_eff  <- sapply(graphs, function(x) mean(brainGraph::efficiency(x)))
global_eff <- sapply(graphs, function(x) brainGraph::efficiency(x, "global"))

# Get statistics for observed networks
obs_df <- GetNetStats(group_fit$networks, group_fit$formula, "model")
colnames(obs_df) <- stat_labels
obs_df <- obs_df %>%
  mutate(n=1:200, Group=rep(c("Young","Old"), each=100)) %>%
  melt(measure.vars=stat_labels, variable.name="Stat")

obs_df <- rbind(obs_df,
                data.frame(n=1:200, Group=rep(c("Young","Old"), each=100), 
                           Stat = "Local efficiency", value = local_eff),
                data.frame(n=1:200, Group=rep(c("Young","Old"), each=100), 
                           Stat = "Global efficiency", value = global_eff))

sim_df <- obs_df[0, ]
young_nets <- group_fit$networks[1:100]
for (i in 1:n_sample){
  y         <- young_nets[[sample(length(young_nets),1)]]
  myformula <- statnet.common::nonsimp_update.formula(group_fit$formula, y ~.,
                                                      from.new = "y")
  
  # Joint fit simul
  this_iter  <- sample(n_iters, 1)
  this_coef <- group_fit$params$muGroup[this_iter,1, ]
  
  net_sim <- simulate(myformula, coef = this_coef,
                     nsim = 1, constraints = constraints,
                     control = control.simulate.formula(MCMC.burnin=aux_iters))
  this_df <- GetNetStats(net_sim, myformula, "model")
  colnames(this_df) <- stat_labels
  this_df <- this_df %>%
    mutate(n=i, Group="Young") %>%
    melt(measure.vars=stat_labels, variable.name="Stat")
  
  this_graph <- igraph::graph_from_adjacency_matrix(as.matrix(net_sim))
  local_eff  <- mean(brainGraph::efficiency(this_graph))
  global_eff <- brainGraph::efficiency(this_graph, "global")
  
  sim_df  <- rbind(sim_df,
                   this_df,
                   data.frame(n=i, Group="Young", 
                              Stat = "Local efficiency", value = local_eff),
                   data.frame(n=i, Group="Young", 
                              Stat = "Global efficiency", value = global_eff))
}

old_nets <- group_fit$networks[101:200]
for (i in 1:n_sample){
  y         <- old_nets[[sample(length(old_nets),1)]]
  myformula <- statnet.common::nonsimp_update.formula(group_fit$formula, y ~.,
                                                      from.new = "y")
  
  # Joint fit simul
  this_iter  <- sample(n_iters, 1)
  this_coef <- group_fit$params$muGroup[this_iter,2, ]

  net_sim <- simulate(myformula, coef = this_coef,
                     nsim = 1, constraints = constraints,
                     control = control.simulate.formula(MCMC.burnin=aux_iters))
  this_df <- GetNetStats(net_sim, myformula, "model")
  colnames(this_df) <- stat_labels
  this_df <- this_df %>%
    mutate(n=i, Group="Old") %>%
    melt(measure.vars=stat_labels, variable.name="Stat")
  
  this_graph <- igraph::graph_from_adjacency_matrix(as.matrix(net_sim))
  local_eff  <- mean(brainGraph::efficiency(this_graph))
  global_eff <- brainGraph::efficiency(this_graph, "global")
  
  sim_df  <- rbind(sim_df,
                   this_df,
                   data.frame(n=i, Group="Old", 
                              Stat = "Local efficiency", value = local_eff),
                   data.frame(n=i, Group="Old", 
                              Stat = "Global efficiency", value = global_eff))
}

obs_df$Observed  <- obs_df$Group
sim_df$Simulated <- sim_df$Group

ggplot(filter(obs_df, Stat %in% c("Local efficiency", "Global efficiency")), 
             aes(Group, value, color = Observed)) + 
  geom_boxplot(width = 0.3, position = position_nudge(x = -0.1, y = 0)) + 
  geom_flat_violin(data = filter(sim_df, Stat %in% c("Local efficiency", "Global efficiency")),
                   colour=NA, alpha = 0.5, width = 0.6, 
                   aes(Group, value, fill=Simulated), inherit.aes = F,
                   position = position_nudge(x = 0.1, y = 0)) +
  facet_wrap("Stat", scales = "free", nrow = 1) +
  scale_fill_discrete(limits = c("Young", "Old"), name = "Posterior predictive") +
  scale_color_discrete(limits = c("Young", "Old"))

ggsave(file.path("Figures", thresh, "gofEfficiency.pdf"), 
       width = 5, height = 2.5)
