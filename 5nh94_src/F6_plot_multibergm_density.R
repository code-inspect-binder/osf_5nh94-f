library(ggplot2)
library(dplyr)
library(cowplot)
theme_set(theme_cowplot(font_size=12))

thresh <- "constK3" # vary between constK3, constK5, propK3, propK5
constraints <- if (substr(thresh,1,1) == "p")  ~edges else ~.

group_fit <- readRDS(file.path("Output", thresh, "combined.RDS"))

# Remove burnIn iterations and apply thinning (default: no thinning)
burn_in        <- 1000
thin          <- 1
post_iters     <- seq(burn_in + 1, group_fit$mainIters, thin)
group_fit$params <- lapply(group_fit$params,
                          function(x) abind::asub(x, post_iters, 1))

stat_labels <- c("edges","gwesp","gwnsp")
if (substr(thresh,1,1) == 'p'){
  stat_labels <- stat_labels[2:3]
}

output <- group_fit$params$muGroup
var_names <- c("Iteration", "Group", "Stat")
out_df <- reshape2::melt(output, varnames = var_names, value.name = "Posterior")
out_df$Group <- factor(out_df$Group, labels = c("Young", "Old"))
out_df$Stat <- factor(out_df$Stat, labels = stat_labels)

ggplot(out_df, aes(x = Posterior)) +
  geom_density(aes(fill = Group), alpha = 0.5, colour = NA) +
  facet_wrap("Stat", ncol = 1, scales = "free") +
  ylab("Density")

ggsave(file.path("Figures", thresh, "groupDiff.pdf"), 
       width = 6.5, height = 3)
