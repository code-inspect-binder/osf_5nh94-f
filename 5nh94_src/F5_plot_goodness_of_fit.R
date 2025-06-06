library(multibergm)
library(ggplot2)
library(dplyr)
library(cowplot)
theme_set(theme_cowplot(font_size=12))

thresh <- "constK3" # vary between constK3, constK5, propK3, propK5
set.seed(42)
constraints <- if (substr(thresh,1,1) == "p")  ~edges else ~.

group_fit <- readRDS(file.path("Output", thresh, "young.RDS"))

# Remove burnIn iterations and apply thinning (default: no thinning)
burn_in       <- 1000
thin          <- 1
post_iters    <- seq(burn_in + 1, group_fit$mainIters, thin)
group_fit$params <- lapply(group_fit$params,
                          function(x) abind::asub(x, post_iters, 1))
n_stats <- length(group_fit$control$model$coef.names)
n_iters <- dim(group_fit$params$muPop)[1]

n_nets    <- length(group_fit$networks)
n_sample  <- 1000
aux_iters <- 100000

# Get statistics for observed networks
obs_df <- GetNetStats(group_fit$networks, group_fit$formula, "gof") %>%
  group_by(Stat) %>%
  mutate(nMax = max(n[Value>0 & is.finite(n)])) %>%
  filter(n <= nMax + 1)

sim_df <- obs_df[0, ]
for (n in 1:n_sample){
  y         <- group_fit$networks[[sample(length(group_fit$networks),1)]]
  myformula <- statnet.common::nonsimp_update.formula(group_fit$formula, y ~.,
                                                      from.new = "y")
  indiv     <- get.network.attribute(y, "individual")
  
  this_iter  <- sample(n_iters, 1)
  this_coef <- group_fit$params$muPop[this_iter, ]
  this_coef <- this_coef + group_fit$params$theta[this_iter,indiv, ]
  net_sim <- simulate(myformula, coef = this_coef,
                     nsim = 1, constraints = constraints,
                     control = control.simulate.formula(MCMC.burnin=aux_iters))
  this_df <- GetNetStats(net_sim, myformula, "gof")
  
  sim_df  <- rbind(sim_df,this_df)
}

sim_df2 <- sim_df %>%
  group_by(Stat, n) %>%
  summarise(Group = mean(Value),
            Lower = quantile(Value, 0.025),
            Upper = quantile(Value, 0.975)) %>%
  mutate(nMax = max(n[Upper>0 & is.finite(n)])) %>%
  filter(n <= nMax + 1)

ggplot(obs_df, aes(x=n, y = Value, group = n)) + 
  geom_boxplot() + 
  geom_line(data = sim_df2, aes(x=n, y = Group), 
            inherit.aes = FALSE) +
  geom_ribbon(data = sim_df2, 
              aes(ymin = Lower, ymax = Upper, x=n), 
              inherit.aes=F,alpha = 0.4) +
  facet_wrap("Stat", scales = "free", ncol = 1) + 
  scale_linetype_discrete(name = NULL, labels = "Observed data") +
  xlab("k") +
  ylab("Probability") + 
  theme(legend.position="top")

ggsave(file.path("Figures", thresh, "gofGroup.pdf"), 
       width = 6.5, height = 6)
