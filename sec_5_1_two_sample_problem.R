# 5. Simulation studies based on real data examples
# 5.1 Two sample problem

hipo <- "H0" # H0, H1
delta <- 1

# for implementation of test
library(fda.usc)
# for Canadian weather data
library(fda)
# for random number generators
library(mvtnorm)

# temperature data
data_set_t <- t(CanadianWeather$dailyAv[,, "Temperature.C"])

# vector of sample sizes
n_vec <- delta * c(15, 15, 5)

# vector of group labels
gr_label_0 <- rep(c(1, 2, 3), c(15, 15, 5))
gr_label <- rep(c(1, 2, 3), n_vec)

nr <- 1000
results <- matrix(NA, nr, 3)
time_1 <- Sys.time()
for (i_nr in seq_len(nr)) {
  if (hipo == "H0") {
    x <- rmvnorm(sum(n_vec), colMeans(data_set_t), cov(data_set_t))
  } else if (hipo == "H1") {
    x_1 <- rmvnorm(n_vec[1], colMeans(data_set_t[gr_label_0 == 1, ]), cov(data_set_t))
    x_2 <- rmvnorm(n_vec[2], colMeans(data_set_t[gr_label_0 == 2, ]), cov(data_set_t))
    x_3 <- rmvnorm(n_vec[3], colMeans(data_set_t[gr_label_0 == 3, ]), cov(data_set_t))
    x <- rbind(x_1, x_2, x_3)
  }
  # matplot(t(x), type = "l", lty = 1, col = rep(1:3, n_vec))
  ec <- x[gr_label == 1, ]
  wc <- x[gr_label == 2, ]
  nc <- x[gr_label == 3, ]
  
  # two sample test for Eastern-Western
  res_1 <- fmean.test.fdata(ec, wc, method = "X2", npc = -0.95)
  # two sample test for Eastern-Northern
  res_2 <- fmean.test.fdata(ec, nc, method = "X2", npc = -0.95)
  # two sample test for Western-Northern
  res_3 <- fmean.test.fdata(wc, nc, method = "X2", npc = -0.95)
  
  results[i_nr, ] <- c(res_1$pvalue[1], res_2$pvalue[1], res_3$pvalue[1])
  print(i_nr)
}
(results_final <- 100 * colMeans(results <= 0.05))
(time_total <- Sys.time() - time_1)

save(results_final, time_total,
     file = paste(paste("z_sim_rde_2", hipo, "d", delta, sep = "_"), "RData", sep = "."))
