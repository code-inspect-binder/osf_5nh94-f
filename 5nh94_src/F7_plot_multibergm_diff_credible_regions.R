library(ggplot2)
library(dplyr)
library(cowplot)
theme_set(theme_cowplot(font_size=12))

thresh <- "constK3" # Vary between constK3 and constK5
constraints <- if (substr(thresh,1,1) == "p")  ~edges else ~.

group_fit <- readRDS(file.path("Output", thresh, "combined.RDS"))

# Remove burnIn iterations and apply thinning (default: no thinning)
burn_in        <- 1000
thin          <- 1
post_iters     <- seq(burn_in + 1, group_fit$mainIters, thin)
group_fit$params <- lapply(group_fit$params,
                          function(x) abind::asub(x, post_iters, 1))

stat_labels <- c("edges","gwesp","gwnsp")

group_posterior <- group_fit$params$muGroup[ ,1, ] - group_fit$params$muGroup[ ,2, ]
colnames(group_posterior) <- stat_labels

# Get single mean fit data
n_terms      <- length(attr(terms(group_fit$formula), "term.labels"))
n_nets <- length(group_fit$networks)/2
singles_mean_all <- array(0, c(15000, n_nets, n_terms))
for (n in 1:n_nets){
  path      <- paste0("Output/", thresh, "/Singles/bergmFit_",sprintf("%03d",n),".RDS")
  single_fit <- readRDS(path)
  singles_mean_all[ ,n, ] <- single_fit$Theta
  
  path      <- paste0("Output/", thresh, "/Singles/bergmFit_",sprintf("%03d",n+100),".RDS")
  single_fit <- readRDS(path)
  singles_mean_all[ ,n, ] <- singles_mean_all[ ,n, ] - single_fit$Theta
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

all_df <- this_df
all_df <- bind_rows(this_df, single_df)
mean_df <- all_df %>%
  group_by(Type) %>%
  summarise_all(mean)

p1 <- ggplot(all_df, aes(edges, gwesp, linetype=Type)) +
  stat_ellipse(type="norm") +
  geom_point(data=mean_df, aes(shape=Type)) +
  geom_vline(xintercept=0, size=0.2) +
  geom_hline(yintercept=0, size=0.2) +
  scale_linetype_discrete(name = NULL, limits = c("multi-BERGM", "mean-BERGM")) +
  scale_shape_manual(name = NULL, limits = c("multi-BERGM", "mean-BERGM"), values = c(16,4))

# gwesp-gwnsp
this_df <- as_tibble(group_posterior[ ,-1])
this_df$Type <- "multi-BERGM"

single_df <- as_tibble(single_sample[ ,-1])
single_df$Type <- "mean-BERGM"

all_df <- this_df
all_df <- bind_rows(this_df, single_df)
mean_df <- all_df %>%
  group_by(Type) %>%
  summarise_all(mean)

p2 <- ggplot(all_df, aes(gwesp, gwnsp, linetype=Type)) +
  stat_ellipse(type="norm") +
  geom_point(data=mean_df, aes(shape=Type)) +
  geom_vline(xintercept=0, size=0.2) +
  geom_hline(yintercept=0, size=0.2) +
  scale_linetype_discrete(name = NULL, limits = c("multi-BERGM", "mean-BERGM")) +
  scale_shape_manual(name = NULL, limits = c("multi-BERGM", "mean-BERGM"), values = c(16,4))

# gwnsp-edges
this_df <- as_tibble(group_posterior[ ,-2])
this_df$Type <- "multi-BERGM"

single_df <- as_tibble(single_sample[ ,-2])
single_df$Type <- "mean-BERGM"

all_df <- this_df
all_df <- bind_rows(this_df, single_df)
mean_df <- all_df %>%
  group_by(Type) %>%
  summarise_all(mean)

p3 <- ggplot(all_df, aes(gwnsp, edges, linetype=Type)) +
  stat_ellipse(type="norm") +
  geom_point(data=mean_df, aes(shape=Type)) +
  geom_vline(xintercept=0, size=0.2) +
  geom_hline(yintercept=0, size=0.2) +
  scale_linetype_discrete(name = NULL, limits = c("multi-BERGM", "mean-BERGM")) +
  scale_shape_manual(name = NULL, limits = c("multi-BERGM", "mean-BERGM"), values = c(16,4))

ggpubr::ggarrange(p1,p2,p3, ncol=3, common.legend = TRUE) 
ggsave(file.path("Figures",thresh, "groupDiffRegion.pdf"), 
       width = 6.5, height = 2.5)

# proportional threshold

thresh <- "propK3" # Vary between propK3 and propK5
group_fit <- readRDS(file.path("Output", thresh, "combined.RDS"))

# Remove burnIn iterations and apply thinning (default: no thinning)
burn_in        <- 1000
thin          <- 1
post_iters     <- seq(burn_in + 1, group_fit$mainIters, thin)
group_fit$params <- lapply(group_fit$params,
                          function(x) abind::asub(x, post_iters, 1))

stat_labels <- c("gwesp","gwnsp")

group_posterior <- group_fit$params$muGroup[ ,1, ] - group_fit$params$muGroup[ ,2, ]
colnames(group_posterior) <- stat_labels

# Get single mean fit data
n_terms      <- length(attr(terms(group_fit$formula), "term.labels"))
n_nets <- length(group_fit$networks)/2
singles_mean_all <- array(0, c(15000, n_nets, n_terms))
for (n in 1:n_nets){
  path      <- paste0("Output/", thresh, "/Singles/bergmFit_",sprintf("%03d",n),".RDS")
  single_fit <- readRDS(path)
  
  singles_mean_all[ ,n, ] <- single_fit$Theta
  singles_mean <- singles_mean + colMeans(single_fit$Theta)
  
  path      <- paste0("Output/", thresh, "/Singles/bergmFit_",sprintf("%03d",n+100),".RDS")
  single_fit <- readRDS(path)
  
  singles_mean_all[ ,n, ] <- singles_mean_all[ ,n, ] - single_fit$Theta
  singles_mean <- singles_mean - colMeans(single_fit$Theta)
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
  geom_vline(xintercept=0, size=0.2) +
  geom_hline(yintercept=0, size=0.2) +
  scale_linetype_discrete(name = NULL, limits = c("multi-BERGM", "mean-BERGM")) +
  scale_shape_manual(name = NULL, limits = c("multi-BERGM", "mean-BERGM"), values = c(16,4))

ggpubr::ggarrange(p1, ncol=1, legend = FALSE)

ggsave(file.path("Figures", thresh, "groupDiffRegion.pdf"), 
       width = 3, height = 2.5)
