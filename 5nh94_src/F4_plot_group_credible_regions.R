library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot(font_size=12))

set.seed(42)

thresh <- "constK3" # vary between constK3 and constK5

# Load in data
output <- readRDS(file.path("Output", thresh, "young.RDS"))

# Remove burnIn iterations and apply thinning (default: no thinning)
burn_in       <- 1000
thin          <- 1
post_iters    <- seq(burn_in + 1, output$mainIters, thin)
output$params <- lapply(output$params,
                         function(x) abind::asub(x, post_iters, 1))

stat_labels <- c("edges","gwesp","gwnsp")

group_posterior <- output$params$muPop
joint_mean      <- colMeans(group_posterior)

colnames(group_posterior) <- stat_labels

n_terms      <- length(attr(terms(output$formula), "term.labels"))
n_nets <- length(output$networks)
singles_mean_all <- array(0, c(15000, n_nets, n_terms))
for (n in 1:n_nets){
  path       <- paste0("Output/", thresh, "/Singles/bergmFit_",sprintf("%03d",n),".RDS")
  single_fit <- readRDS(path)
  
  singles_mean_all[ ,n, ] <- single_fit$Theta
}
singles_mu <- apply(singles_mean_all, c(1,3), mean)
singles_mean <- colMeans(singles_mu)
singles_cov <- cov(singles_mu) + 
  matrix(rowMeans(apply(singles_mean_all, 1, var)/n_nets), n_terms,n_terms)

single_sample <- mvtnorm::rmvnorm(30000, singles_mean, singles_cov)
colnames(single_sample) <- stat_labels

# edges-gwesp
this_df <- as_tibble(group_posterior[ ,-3])
this_df$Type <- "multi-BERGM"

single_df <- as_tibble(single_sample[ ,-3])
single_df$Type <- "mean-BERGM"

all_df <- bind_rows(this_df, single_df)
mean_df <- all_df %>%
  group_by(Type) %>%
  summarise_all(mean)

p1 <- ggplot(all_df, aes(edges, gwesp, linetype=Type)) +
  stat_ellipse(type="norm") +
  geom_point(data=mean_df, aes(shape=Type)) +
  scale_linetype_discrete(name = NULL, limits = c("multi-BERGM", "mean-BERGM")) +
  scale_shape_manual(name = NULL, limits = c("multi-BERGM", "mean-BERGM"), values = c(16,4))

# gwesp-gwnsp
this_df <- as_tibble(group_posterior[ ,-1])
this_df$Type <- "multi-BERGM"

single_df <- as_tibble(single_sample[ ,-1])
single_df$Type <- "mean-BERGM"

all_df <- bind_rows(this_df, single_df)
mean_df <- all_df %>%
  group_by(Type) %>%
  summarise_all(mean)

p2 <- ggplot(bind_rows(this_df, single_df), aes(gwesp, gwnsp, linetype=Type)) +
  stat_ellipse(type="norm") +
  geom_point(data=mean_df, aes(shape=Type)) +
  scale_linetype_discrete(name = NULL, limits = c("multi-BERGM", "mean-BERGM")) +
  scale_shape_manual(name = NULL, limits = c("multi-BERGM", "mean-BERGM"), values = c(16,4))

# gwnsp-edges
this_df <- as_tibble(group_posterior[ ,-2])
this_df$Type <- "multi-BERGM"

single_df <- as_tibble(single_sample[ ,-2])
single_df$Type <- "mean-BERGM"

all_df <- bind_rows(this_df, single_df)
mean_df <- all_df %>%
  group_by(Type) %>%
  summarise_all(mean)

p3 <- ggplot(bind_rows(this_df, single_df), aes(gwnsp, edges, linetype=Type)) +
  stat_ellipse(type="norm") +
  geom_point(data=mean_df, aes(shape=Type)) +
  scale_linetype_discrete(name = NULL, limits = c("multi-BERGM", "mean-BERGM")) +
  scale_shape_manual(name = NULL, limits = c("multi-BERGM", "mean-BERGM"), values = c(16,4))

ggpubr::ggarrange(p1,p2,p3, ncol=3, common.legend = TRUE) 
ggsave(file.path("Figures", thresh, "credibleRegions.pdf"), 
       width = 6.5, height = 2.5)

#####
# proportional threshold plots
set.seed(42)
thresh <- "propK3" # Vary between propK3 and propK5

# Load in data
output <- output <- readRDS(file.path("Output", thresh, "young.RDS"))

# Remove burnIn iterations and apply thinning (default: no thinning)
burn_in       <- 1000
thin          <- 1
post_iters    <- seq(burn_in + 1, output$mainIters, thin)
output$params <- lapply(output$params,
                        function(x) abind::asub(x, post_iters, 1))

stat_labels <- c("gwesp","gwnsp")

group_posterior <- output$params$muPop
joint_mean      <- colMeans(group_posterior)

colnames(group_posterior) <- stat_labels

n_terms      <- length(attr(terms(output$formula), "term.labels"))
n_nets <- length(output$networks)
singles_mean_all <- array(0, c(15000, n_nets, n_terms))
for (n in 1:n_nets){
  path      <- paste0("Output/", thresh, "/Singles/bergmFit_",sprintf("%03d",n),".RDS")
  single_fit <- readRDS(path)
  
  singles_mean_all[ ,n, ] <- single_fit$Theta
}
singles_mu <- apply(singles_mean_all, c(1,3), mean)
singles_mean <- colMeans(singles_mu)
singles_cov <- cov(singles_mu) + 
  matrix(rowMeans(apply(singles_mean_all, 1, var)/n_nets), n_terms,n_terms)

single_sample <- mvtnorm::rmvnorm(30000, singles_mean, singles_cov)
colnames(single_sample) <- stat_labels

# gwesp-gwnsp
this_df <- as_tibble(group_posterior)
this_df$Type <- "multi-BERGM"

single_df <- as_tibble(single_sample)
single_df$Type <- "mean-BERGM"

all_df <- bind_rows(this_df, single_df)
mean_df <- all_df %>%
  group_by(Type) %>%
  summarise_all(mean)

p1 <- ggplot(all_df, aes(gwesp, gwnsp, linetype=Type)) +
  stat_ellipse(type="norm") +
  geom_point(data=mean_df, aes(shape=Type)) +
  scale_linetype_discrete(name = NULL, limits = c("multi-BERGM", "mean-BERGM")) +
  scale_shape_manual(name = NULL, limits = c("multi-BERGM", "mean-BERGM"), values = c(16,4))

ggpubr::ggarrange(p1, ncol=1, legend="none")

ggsave(file.path("Figures",thresh, "credibleRegions.pdf"), 
       width = 3, height = 2.5)
