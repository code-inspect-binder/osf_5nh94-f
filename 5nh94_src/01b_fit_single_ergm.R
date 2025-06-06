### IMPORTANT: Install forked version of Bergm to deal with constraints
# devtools::install_github("brieuclehmann/Bergm")
library(Bergm)

param  <- commandArgs(trailingOnly = TRUE)
thresh <- param[1]
i  <- as.integer(param[2])
set.seed(42 + i) # for reproducibility

if (!(thresh %in% c("constK3", "constK5", "propK3", "propK5"))) {
  stop("Incorrect threshold specified")
}

# Get network
y <- readRDS(paste0("Data/n200_", thresh, ".RDS"))[[i]]

# Set up ERGM formula and constraints
if (substr(thresh, 1, 4) == 'prop') {
  ergm_formula <- y ~ gwesp(0.75, fixed=TRUE) + gwnsp(0.75, fixed=TRUE)
  constraints <- ~edges
} else {
  ergm_formula <- y ~ edges + gwesp(0.75, fixed=TRUE) + gwnsp(0.75, fixed=TRUE)
  constraints <- ~.
}

ergm_fit <- ergm(ergm_formula, constraints = constraints)

f <- paste0("Output/Singles/", thresh, "/ergmFit_", sprintf("%03d",i), ".RDS")
dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
saveRDS(bergm_fit, f)
