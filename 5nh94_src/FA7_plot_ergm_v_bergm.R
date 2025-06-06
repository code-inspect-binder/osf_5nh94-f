set.seed(42) # for reproducibility
library(ergm)
library(tibble)

thresh <- "constK3"
net_ind <- 1:200

# Set up ERGM formula and constraints
if (substr(thresh,1,1) == 'p') {
  ergm_formula <- nets ~ gwesp(0.75, fixed=TRUE) + gwnsp(0.75, fixed=TRUE)
  constraints <- ~edges
} else {
  ergm_formula <- nets ~ edges + gwesp(0.75, fixed=TRUE) + gwnsp(0.75, fixed=TRUE)
  constraints <- ~.
}
n_terms <- length(attr(terms(ergm_formula),'term.labels'))

bergm_files <- paste0("Output/", thresh, "/Singles/bergmFit_",sprintf("%03d",net_ind),".RDS")
all_theta <- sapply(bergm_files, function(x) readRDS(x)$Theta, 
                    simplify = "array", USE.NAMES = FALSE)
bergm_theta <- apply(all_theta, c(2,3), mean)

ergm_files <- list.files(paste0("Output/", thresh, "/Singles/"), '^ergmFit', 
                         full.names = TRUE)

# Manual inspection of mcmc.diagnostics() plots:
ind_manual <- c(20, 26, 31, 36, 37, 39, 42, 58, 63, 67, 86,
                111, 112, 140, 149, 155, 162, 170, 174)
did_not_converge <- sapply(ind_manual, 
                           function(x) as.integer(substr(ergm_files[x],41,43)))
  
ergm_theta <- matrix(NA, 3, 200)
msngno <- vector("integer")        
for (i in 1:200) {
  f <- paste0("Output/", thresh, "/Singles/ergmFit_", sprintf("%03d", i), ".RDS")
  if (!file.exists(f)) {
    msngno <- c(msngno, i)
  } else {
    ergm_theta[ ,i] <- coef(readRDS(f))
  }
}

ergm_theta[ ,did_not_converge] <- NA

par(mfrow = c(1,3))
plot(ergm_theta[1, ], bergm_theta[1, ], col = c(rep(1,100), rep(2,100)),
     xlab = "edges (ergm)", ylab = "edges (bergm)")
abline(0, 1)
plot(ergm_theta[2, ], bergm_theta[2, ], col = c(rep(1,100), rep(2,100)),
     xlab = "gwesp (ergm)", ylab = "gwesp (bergm)")
abline(0, 1)
plot(ergm_theta[3, ], bergm_theta[3, ], col = c(rep(1,100), rep(2,100)),
     xlab = "gwnsp (ergm)", ylab = "gwnsp (bergm)")
abline(0, 1)

set.seed(1) # set seed for reproducibility
n_perm <- 100000
null_dist <- matrix(NA, n_perm, 3)
p_value <- rep(NA, 3)
t_score <- rep(NA, 3)
obs_diff <- rep(NA, 3)

for (j in 1:3){
  for (i in 1:n_perm) {
    shuffled_data <- sample(bergm_theta[j, ])
    shuffled_young <- shuffled_data[1:100] 
    shuffled_old <- shuffled_data[101:200] 
    null_dist[i, j] <- mean(shuffled_young) - mean(shuffled_old)
  }
  obs_diff[j] <- mean(bergm_theta[j,1:100]) - mean(bergm_theta[j,101:200])
  sd_tot <- sqrt(var(bergm_theta[j,1:100]) + var(bergm_theta[j,101:200]))
  t_score[j] <- obs_diff[j]/sd_tot
  p_value[j] <- mean(abs(null_dist[, j]) >= abs(obs_diff[j]))
}

ttest_bergm_df <- tibble(param = c("edges", "gwesp", "gwnsp"),
                         diff = obs_diff,
                         p_value = p_value)

set.seed(1) # set seed for reproducibility
n_perm <- 100000
null_dist <- matrix(NA, n_perm, 3)
p_value <- rep(NA, 3)
t_score <- rep(NA, 3)
obs_diff <- rep(NA, 3)

ind_young <- which(!is.na(ergm_theta[1,1:100]))
ind_old <- 100 + which(!is.na(ergm_theta[1,101:200]))

for (j in 1:3){
  for (i in 1:n_perm) {
    shuffled_data <- sample(ergm_theta[j, c(ind_young, ind_old)])
    shuffled_young <- shuffled_data[1:length(ind_young)] 
    shuffled_old <- shuffled_data[-(1:length(ind_young))] 
    null_dist[i, j] <- mean(shuffled_young) - mean(shuffled_old)
  }
  obs_diff[j] <- mean(ergm_theta[j,ind_young]) - mean(ergm_theta[j,ind_old])
  sd_tot <- sqrt(var(ergm_theta[j,ind_young]) + var(ergm_theta[j,ind_old]))
  t_score[j] <- obs_diff[j]/sd_tot
  p_value[j] <- mean(abs(null_dist[, j]) >= abs(obs_diff[j]))
}

ttest_ergm_df <- tibble(param = c("edges", "gwesp", "gwnsp"),
                        diff = obs_diff,
                        p_value = p_value)
