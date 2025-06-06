nets <- readRDS("Data/n200_constK3.RDS")

model_stats <- sapply(nets, function(x) summary(x ~ edges + gwesp(0.75, fixed=TRUE) + gwnsp(0.75, fixed=TRUE)))
glob_eff <- sapply(nets, function(x) brainGraph::efficiency(igraph::graph_from_adjacency_matrix(as.matrix(x)), "global"))
loc_eff <- sapply(nets, function(x) mean(brainGraph::efficiency(igraph::graph_from_adjacency_matrix(as.matrix(x)), "local")))

all_stats <- rbind(model_stats, loc_eff, glob_eff)
all_stats <- t(all_stats)

xtable(cor(all_stats[ ,c(2,3,4,5)], method = "spearman"))

xtable(ppcor::pcor(all_stats[ ,c(2,3,4)], method = "spearman")$estimate)
xtable(ppcor::pcor(all_stats[ ,c(2,3,4)], method = "spearman")$p.value)

xtable(ppcor::pcor(all_stats[ ,c(2,3,5)], method = "spearman")$estimate)
xtable(ppcor::pcor(all_stats[ ,c(2,3,5)], method = "spearman")$p.value)
