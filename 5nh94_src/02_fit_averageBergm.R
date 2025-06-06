library(Bergm)
set.seed(42)

param  <- commandArgs(trailingOnly=TRUE)
thresh <- param[1]
type   <- param[2]

if (!(thresh %in% c("constK3", "constK5", "propK3", "propK5"))) {
  stop("Incorrect threshold specified")
}

if (!(type %in% c("mean", "median"))) {
  stop("Incorrect type specified")
}

y <- readRDS(paste0("Data/n200_young_", param, ".RDS"))
ergm_formula <- y ~ edges + gwesp(0.75, fixed = TRUE) + gwnsp(0.75, fixed = TRUE)

# Vary this proposal parameter to get average acceptance rate of around 0.3
sigma_proposal <- 0.003
bergm_fit <- bergm(ergm_formula, main.iters = 5000, nchains = 6,
                   aux.iters = 20000, burn.in = 1000, 
                   V.proposal = sigma_proposal)

out_file <- file.path("Output", thresh, paste0(type, "_young.RDS"))
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(bergm_fit, out_file)
