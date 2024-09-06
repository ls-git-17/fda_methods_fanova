# Supplement
# 1. Simulation studies based on real data examples
# 1.4 Multivariate analysis of variance

hipo <- "H0" # H0, H1
delta <- 1

# for FMANOVA test by Gorecki and Smaga (2017)
library(fdANOVA)
# for Canadian weather data
library(fda)
# for random number generators
library(mvtnorm)

# temperature and precipitation data
data_set_t <- t(CanadianWeather$dailyAv[,, "Temperature.C"])
data_set_p <- t(CanadianWeather$dailyAv[,, "Precipitation.mm"])
data_set_t_p_b <- cbind(data_set_t, data_set_p)
# data for both features
data_set_t_p <- vector("list", 2)
data_set_t_p[[1]] <- t(data_set_t)
data_set_t_p[[2]] <- t(data_set_p)

# vector of sample sizes
n_vec <- delta * c(15, 15, 5)

# vector of group labels
gr_label_0 <- rep(c(1, 2, 3), c(15, 15, 5))
gr_label <- rep(c(1, 2, 3), n_vec)

nr <- 1000
results <- matrix(NA, nr, 28)
time_1 <- Sys.time()
for (i_nr in seq_len(nr)) {
  print(i_nr)
  if (hipo == "H0") {
    temp_h0 <- rmvnorm(sum(n_vec), colMeans(data_set_t_p_b), cov(data_set_t_p_b))
    data_set_t_p_sim <- vector("list", 2)
    data_set_t_p_sim[[1]] <- t(temp_h0[, 1:365])
    data_set_t_p_sim[[2]] <- t(temp_h0[, 366:730])
  } else if (hipo == "H1") {
    temp_h1_1 <- rmvnorm(n_vec[1], colMeans(data_set_t_p_b[gr_label_0 == 1, ]), cov(data_set_t_p_b))
    temp_h1_2 <- rmvnorm(n_vec[2], colMeans(data_set_t_p_b[gr_label_0 == 2, ]), cov(data_set_t_p_b))
    temp_h1_3 <- rmvnorm(n_vec[3], colMeans(data_set_t_p_b[gr_label_0 == 3, ]), cov(data_set_t_p_b))
    data_set_t_p_sim <- vector("list", 2)
    data_set_t_p_sim[[1]] <- t(rbind(temp_h1_1[, 1:365], temp_h1_2[, 1:365], temp_h1_3[, 1:365]))
    data_set_t_p_sim[[2]] <- t(rbind(temp_h1_1[, 366:730], temp_h1_2[, 366:730], temp_h1_3[, 366:730]))
  }
  # matplot(t(temp_h0), type = "l", lty = 1, col = rep(1:3, n_vec))
  # matplot(data_set_t_p_sim[[1]], type = "l", lty = 1, col = rep(1:3, n_vec))
  # matplot(data_set_t_p_sim[[2]], type = "l", lty = 1, col = rep(1:3, n_vec))
  # matplot(t(temp_h1_1), type = "l", lty = 1, col = rep(1:3, n_vec))
  # matplot(t(temp_h1_2), type = "l", lty = 1, col = rep(1:3, n_vec))
  # matplot(t(temp_h1_3), type = "l", lty = 1, col = rep(1:3, n_vec))
  
  # the tests based on a basis function representation
  fmanova_ptbfr <- fmanova.ptbfr(data_set_t_p_sim, gr_label, maxK = 31)
  
  # the tests based on random projections with the Gaussian white noise
  fmanova_trp_1 <- fmanova.trp(data_set_t_p_sim, gr_label, k = c(10, 20, 30))
  
  # the tests based on random projections with the Brownian motion
  fmanova_trp_2 <- fmanova.trp(data_set_t_p_sim, gr_label, k = c(10, 20, 30), 
                               projection = "BM")
  
  results[i_nr, ] <- c(fmanova_ptbfr$pvalueW, fmanova_ptbfr$pvalueLH, fmanova_ptbfr$pvalueP, fmanova_ptbfr$pvalueR,
                       fmanova_trp_1$pvalues[1, ], fmanova_trp_1$pvalues[2, ], fmanova_trp_1$pvalues[3, ], fmanova_trp_1$pvalues[4, ],
                       fmanova_trp_2$pvalues[1, ], fmanova_trp_2$pvalues[2, ], fmanova_trp_2$pvalues[3, ], fmanova_trp_2$pvalues[4, ])
}
(results_final <- 100 * colMeans(results <= 0.05))
(time_total <- Sys.time() - time_1)

save(results_final, time_total,
     file = paste(paste("z_sim_rde_5", hipo, "d", delta, sep = "_"), "RData", sep = "."))

# the tests based on a basis function representation
results_final[1:4]
# the tests based on random projections with the Gaussian white noise
# W
results_final[5:7]
# LH
results_final[8:10]
# P
results_final[11:13]
# R
results_final[14:16]
# the tests based on random projections with the Brownian motion
# W
results_final[17:19]
# LH
results_final[20:22]
# P
results_final[23:25]
# R
results_final[26:28]

plot(results_final, type = "l")
text(1:28, results_final, 1:28)
