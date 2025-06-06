library(dplyr)
library(ggplot2)
library(patchwork)
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

net_ind <- c(10, 53, 78)
pl <- list()

for(i in 1:length(net_ind)){
  s <- net_ind[i]
  group_posterior <- output$params$theta[ ,s, ] + output$params$muPop
  
  colnames(group_posterior) <- stat_labels
  
  path       <- paste0("Output/", thresh, "/Singles/bergmFit_",sprintf("%03d",s),".RDS")
  single_fit <- readRDS(path)
  
  single_sample <- single_fit$Theta
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
  
  p1 <- ggplot(bind_rows(this_df, single_df), aes(edges, gwesp, linetype=Type)) +
    stat_ellipse(type="norm") +
    geom_point(data=mean_df, aes(shape=Type)) +
    scale_linetype_discrete(name = NULL, limits = c("multi-BERGM", "mean-BERGM")) +
    scale_shape_manual(name = NULL, limits = c("multi-BERGM", "mean-BERGM"), values = c(16,4)) +
    ggtitle(paste("i =", s))
  
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
  
  #if (i == 2)
  #  pl[[i]] <- ggpubr::ggarrange(p1,p2,p3, nrow=3, common.legend=T)
  #else
    pl[[i]] <- ggpubr::ggarrange(p1,p2,p3, nrow=3, legend=F)
}
ggpubr::ggarrange(pl[[1]],pl[[2]],pl[[3]], ncol=3, legend=F)

ggsave(file.path("Figures", thresh, "indivRegions.pdf"), 
       width = 6, height = 5)

#####
# proportional threshold plots
set.seed(42)
thresh <- "propK3" # vary between propK3 and propK5

# Load in data
output <- readRDS(file.path("Output", thresh, "young.RDS"))

# Remove burnIn iterations and apply thinning (default: no thinning)
burn_in       <- 1000
thin          <- 1
post_iters    <- seq(burnIn + 1, output$mainIters, thin)
output$params <- lapply(output$params,
                        function(x) abind::asub(x, postIters, 1))

stat_labels <- c("gwesp","gwnsp")

net_ind <- c(10, 53, 78)
pl <- list()

for(i in 1:length(net_ind)){
  s <- net_ind[i]
  group_posterior <- output$params$theta[ ,s, ] + output$params$muPop
  
  colnames(group_posterior) <- stat_labels
  
  path       <- paste0("Output/", thresh, "/Singles/bergmFit_",sprintf("%03d",s),".RDS")
  single_fit <- readRDS(path)
  
  single_sample <- single_fit$Theta
  colnames(single_sample) <- stat_labels
  
  # gwesp-gwnsp
  this_df <- as_tibble(group_posterior)
  this_df$Type <- "multi-BERGM"
  
  single_df <- as_tibble(single_sample[ , ])
  single_df$Type <- "mean-BERGM"
  
  all_df <- bind_rows(this_df, single_df)
  mean_df <- all_df %>%
    group_by(Type) %>%
    summarise_all(mean)
  
  p1 <- ggplot(bind_rows(this_df, single_df), aes(gwesp, gwnsp, linetype=Type)) +
    stat_ellipse(type="norm") +
    geom_point(data=mean_df, aes(shape=Type)) +
    scale_linetype_discrete(name = NULL, limits = c("multi-BERGM", "mean-BERGM")) +
    scale_shape_manual(name = NULL, limits = c("multi-BERGM", "mean-BERGM"), values = c(16,4)) +
    ggtitle(paste("i =", s))
  
  #if (i == 2)
  #  pl[[i]] <- ggpubr::ggarrange(p1,p2,p3, nrow=3, common.legend=T)
  #else
  pl[[i]] <- ggpubr::ggarrange(p1, nrow=1, legend=F)
}
ggpubr::ggarrange(pl[[1]],pl[[2]],pl[[3]], ncol=3, legend=F)
ggsave(file.path("Figures", thresh, "indivRegions.pdf"), 
       width = 6, height = 2)
