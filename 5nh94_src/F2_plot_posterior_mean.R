library(cowplot)
library(ggplot)
library(dplyr)
theme_set(theme_cowplot(font_size=12))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n_nets <- 100
thresh <- "constK3" # vary between constK3, constK5, propK3, propK5

mean_fit <- readRDS(file.path("Output", thresh, "mean_young.RDS"))
mean_post_mean <- colMeans(mean_fit$Theta)

median_fit <- readRDS(file.path("Output", thresh, "median_young.RDS"))
median_post_mean <- colMeans(median_fit$Theta)

group_fit <- readRDS(file.path("Output", thresh, "young.RDS"))
# Remove burnIn iterations and apply thinning (default: no thinning)
burn_in        <- 1000
thin          <- 1
post_iters     <- seq(burn_in + 1, group_fit$mainIters, thin)
group_fit$params <- lapply(group_fit$params,
                          function(x) abind::asub(x, post_iters, 1))
group_post_mean <- colMeans(group_fit$params$muPop)

group_single_mean <- t(apply(group_fit$params$theta, 2, colMeans) + group_post_mean)

n_stats <- length(group_post_mean)
single_post_mean <- matrix(NA, n_nets, n_stats)
for (n in 1:n_nets){
  f <- paste0("Output/", thresh, "/Singles/bergmFit_",sprintf("%03d",n),".RDS")
  singleFit <- readRDS(f)
  
  single_post_mean[n, ] <- colMeans(singleFit$Theta)
}
mean_singles <- colMeans(single_post_mean)

means_df <- tibble(Type=character(), Stat = character(), Mean = double())
singles_df <- tibble(Type=character(), Stat = character(), Mean = double())
for (i in 1:n_stats){
  out_df <- tibble(Mean = c(mean_singles[i], 
                            group_post_mean[i],
                            mean_post_mean[i], 
                            median_post_mean[i]), 
                   Type= c("mean-BERGM", "multi-BERGM", 
                           "mean-GRN BERGM", "median-GRN BERGM"),
                   Stat = group_fit$control$model$coef.names[i])
  
  means_df <- rbind(means_df,out_df)
  
  out_df <- tibble(Mean = single_post_mean[ ,i],
                   Individual = "",
                   Stat = group_fit$control$model$coef.names[i])
  
  singles_df <- rbind(singles_df,out_df)
}

if (thresh %in% c("constK3", "constK5")) {
  
  singles_df <- singles_df %>%
    mutate(Stat = factor(Stat, labels = c("edges","gwesp","gwnsp")))
  
  means_df <- means_df %>%
    mutate(Stat = factor(Stat, labels = c("edges","gwesp","gwnsp")))
  
  ggplot(singles_df, aes(x = Mean, fill = Individual)) +
    geom_histogram(data = filter(singles_df, Stat == "edges"),
                   binwidth = 0.6, alpha = 0.5, position = "identity") +
    geom_histogram(data = filter(singles_df, Stat == "gwesp"),
                   binwidth = 0.2, alpha = 0.5, position = "identity") +
    geom_histogram(data = filter(singles_df, Stat == "gwnsp"),
                   binwidth = 0.1, alpha = 0.5, position = "identity") +
    facet_wrap("Stat", scales = "free_x") +
    geom_vline(data = means_df, 
               aes(xintercept = Mean, color = Type, linetype = Type), 
               size = 1) +
    scale_fill_manual(values = "gray", name = NULL, 
                      labels = "Individual fits") +
    scale_color_manual(values = gg_color_hue(4)[c(1,3,2,4)], limits = c("multi-BERGM", "mean-BERGM", "mean-GRN BERGM", "median-GRN BERGM"),
                       name = "Group-based fits") +
    scale_linetype_discrete(name = "Group-based fits", limits = c("multi-BERGM", "mean-BERGM", "mean-GRN BERGM", "median-GRN BERGM")) +
    scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
    xlab("Posterior mean") + 
    ylab("Count")
  
  ggsave(file.path("Figures",thresh,"posteriorMean.pdf"), 
         width = 6.5, height = 2)
  
} else {
  singles_df <- rbind(c(NA,"","edges"), singles_df)
  singles_df <- singles_df %>%
    mutate(Mean = as.double(Mean)) %>%
    mutate(Stat = factor(Stat, levels = unique(Stat), 
                         labels = c("edges","gwesp","gwnsp")))
  
  means_df <- rbind(c(NA,"mean-GRN BERGM","edges"), means_df)
  means_df <- means_df %>%
    mutate(Mean = as.double(Mean)) %>%
    mutate(Stat = factor(Stat, levels = unique(Stat), 
                         labels = c("edges", "gwesp","gwnsp")))
  
  ggplot(singles_df, aes(x = Mean, fill = Individual)) +
    geom_histogram(data = filter(singles_df, Stat == "gwesp"), boundary = 1,
                   binwidth = 0.2, alpha = 0.5, position = "identity") +
    geom_histogram(data = filter(singles_df, Stat == "gwnsp"), boundary = 0.5,
                   binwidth = 0.1, alpha = 0.5, position = "identity") +
    geom_vline(data = means_df, 
               aes(xintercept = Mean, color = Type, linetype = Type), 
               size = 1.2) +
    facet_wrap("Stat", scales = "free_x") +
    scale_fill_manual(values = "gray", name = NULL, 
                      labels = "Individual fits") +
    scale_color_manual(values = gg_color_hue(4)[c(1,3,2,4)], limits = c("multi-BERGM", "mean-BERGM", "mean-GRN BERGM", "median-GRN BERGM"),
                       name = "Group-based fits") +
    scale_linetype_discrete(name = "Group-based fits", limits = c("multi-BERGM", "mean-BERGM", "mean-GRN BERGM", "median-GRN BERGM")) +
    xlab("Posterior mean") + 
    ylab("Count")
  
  ggsave(file.path("Figures",thresh,"posteriorMean.pdf"), 
         width = 6.5, height = 2)
}
