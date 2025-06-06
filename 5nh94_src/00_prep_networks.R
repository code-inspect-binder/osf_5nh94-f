### Define network construction function
library(network)
constructNetworks <- function(fc, threshold = 0.1, proportion = NULL){
  
  n_nets  <- dim(fc)[1]
  n_nodes <- dim(fc)[2]
  if (is.null(proportion)) {
    n_edges <- floor(threshold * n_nodes * (n_nodes - 1) / 2)
    cutoff <- sort(fc, decreasing = TRUE)[2 * n_edges * n_nets]
  }
  nets <- list()
  for (n in 1:n_nets){
    if (!is.null(proportion)) {
      n_edges <- floor(proportion*n_nodes*(n_nodes-1)/2)
      cutoff <- sort(fc[n, , ], decreasing=TRUE)[2*n_edges]
    }
    this_mat <- fc[n, , ]
    this_net <- network(this_mat >= cutoff, directed = FALSE)
    this_net <- set.network.attribute(this_net, attrname = "individual",
                                      value = n)
    nets <- append(nets, list(this_net))
  }
  nets
}

### Construct and save networks
fc <- readRDS("Data/fc_mats.RDS")

# Get indices of 100 youngest and 100 oldest participants
n_sub     <- dim(fc)[1]
young_ind <- seq(1, 100)
old_ind   <- seq(n_sub - 100 + 1, n_sub)
net_ind   <- c(young_ind, old_ind)

n_nodes <- dim(fc)[2]
avg_node_degree1 <- 3 # ECO thresholding
avg_node_degree2 <- exp(log(n_nodes) / 2.8) # Simpson et al. (2011)

for (avg_node_degree in c(avg_node_degree1, avg_node_degree2)) {
  # Calculate threshold based on average node degree
  threshold <- avg_node_degree / (n_nodes - 1)
  
  nets <- constructNetworks(fc[net_ind, , ], threshold = threshold)
  saveRDS(nets, paste0("Data/n200_constK", round(avg_node_degree), ".RDS"))
  
  nets <- constructNetworks(fc[net_ind, , ], proportion = threshold)
  saveRDS(nets, paste0("Data/n200_propK", round(avg_node_degree), ".RDS"))
  
  # Construct mean-representative networks
  for (avg in c("mean", "median")) {
    avg_mat <- array(apply(fc[young_ind, , ], c(2, 3), avg), 
                     dim = c(1, n_nodes, n_nodes))
    nets    <- constructNetworks(avg_mat, proportion = threshold)
    saveRDS(nets[[1]], paste0("Data/n200_young_", avg, "K", 
                              round(avg_node_degree), ".RDS"))
    
    avg_mat <- array(apply(fc[old_ind, , ], c(2, 3), avg), 
                     dim = c(1, n_nodes, n_nodes))
    nets    <- constructNetworks(avg_mat, proportion = threshold)
    saveRDS(nets[[1]], paste0("Data/n200_old_", avg, "K", 
                              round(avg_node_degree), ".RDS"))
  }
}