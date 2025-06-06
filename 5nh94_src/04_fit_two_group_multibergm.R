# IMPORTANT: install release v0.1 of multibergm
# devtools::install_github("brieuclehmann/multibergm@v0.1")
library(multibergm)

param          <- commandArgs(trailingOnly=TRUE)
thresh         <- param[1]
theta_proposal <- as.double(param[2])
mu_proposal    <- as.double(param[3])
main_iters     <- as.integer(param[4])
n_batches      <- as.double(param[5])

set.seed(42) # for reproducibility

# Get networks
young_ind       <- 1:100
old_ind         <- 101:200
net_ind         <- c(young_ind, old_ind)
nets            <- readRDS(paste0("Data/n200_", thresh, ".RDS"))
nets[young_ind] <- lapply(nets[young_ind], 
                         function(x) set.network.attribute(x, "group", 1))
nets[old_ind]   <- lapply(nets[old_ind], 
                         function(x) set.network.attribute(x, "group", 2))

print(paste("Running threshold",thresh))

# Set up ERGM formula and constraints
if (substr(thresh,1,1) == 'p') {
  ergm_formula <- nets ~ gwesp(0.75, fixed=TRUE) + gwnsp(0.75, fixed=TRUE)
  constraints <- ~edges
} else {
  ergm_formula <- nets ~ edges + gwesp(0.75, fixed=TRUE) + gwnsp(0.75, fixed=TRUE)
  constraints <- ~.
}
n_terms <- length(attr(terms(ergm_formula),'term.labels'))

# Get posterior values from individual fits
all_files <- paste0("Output/Singles/", thresh, "/bergmFit_", 
                   sprintf("%03d",net_ind), ".RDS")
all_theta <- sapply(all_files, function(x) readRDS(x)$Theta, 
                    simplify = "array", USE.NAMES = FALSE)

singles_mean <- t(apply(all_theta, c(2,3), mean))
group_mean   <- rbind(colMeans(singles_mean[young_ind, ]),
                     colMeans(singles_mean[old_ind, ]))
pop_mean     <- colMeans(singles_mean)
output_cov   <- array(NA, c(length(nets), n_terms, n_terms))
for (n in net_ind)
  output_cov[n, , ] <- cov(all_theta[ , ,n])
group_cov    <- apply(output_cov, c(2,3), sum)

# Set up initial values
init          <- list()
init$theta    <- array(sweep(singles_mean,2,pop_mean), c(length(nets),n_terms))
init$muPop    <- pop_mean
init$covTheta <- array(cov(singles_mean), c(1,n_terms,n_terms))

# Set up proposals
proposal <- list(theta = group_cov*theta_proposal, 
                 mu    = cov(singles_mean)*mu_proposal)

# Set up priors
df_cov <- n_terms
prior <- list(theta   = list(scale = df_cov*cov(singles_mean),
                             df = df_cov),
              muGroup = list(scale = df_cov*cov(singles_mean),
                             df = df_cov),
              muPop   = list(mean = pop_mean,
                             cov = group_cov))

control <- control_multibergm(ergm_formula, proposal = proposal, 
                              constraints = constraints,
                              prior = prior, init = init, 
                              nBatches=n_batches, auxIters=20000)

# Fit multibergm
fit <- multibergm(ergm_formula, constraints, main_iters, control)

out_file <- file.path("Output", thresh, "combined.RDS")
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(fit, outFile)
