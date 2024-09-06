# Supplement
# 1. Simulation studies based on real data examples
# 1.3 ANOVA for partially observed functional data

hipo <- "H0" # H0, H1
delta <- 1

# for implementation of test
library(fda.usc)
# for Canadian weather data
library(fda)
# for random number generators
library(mvtnorm)
# for data storing
library(data.table)

# dowload the supplementary code for the paper
# Kraus D. (2019). Inferential procedures for partially observed functional data.
# Journal of Multivariate Analysis 173, 583-603
# from the webpage (available 2023.10.11)
# https://www.sciencedirect.com/science/article/pii/S0047259X18304950#appSD

# from downloaded files load the following scripts:
# for simul.obs() function
source("simul2.fd.R")
# for meantest() function
source("mean.test.missfd.R")

# temperature data
data_set_t <- t(CanadianWeather$dailyAv[,, "Temperature.C"])

# vector of sample sizes
n_vec <- delta * c(15, 15, 5)

# vector of group labels
gr_label_0 <- rep(c(1, 2, 3), c(15, 15, 5))
gr_label <- rep(c(1, 2, 3), n_vec)

# generate partially observed functional data from temperature data
ti <- 1:365

nr <- 1000
results <- matrix(NA, nr, 3)
time_1 <- Sys.time()
for (i_nr in seq_len(nr)) {
  print(i_nr)
  if (hipo == "H0") {
    x <- rmvnorm(sum(n_vec), colMeans(data_set_t), cov(data_set_t))
  } else if (hipo == "H1") {
    x_1 <- rmvnorm(n_vec[1], colMeans(data_set_t[gr_label_0 == 1, ]), cov(data_set_t))
    x_2 <- rmvnorm(n_vec[2], colMeans(data_set_t[gr_label_0 == 2, ]), cov(data_set_t))
    x_3 <- rmvnorm(n_vec[3], colMeans(data_set_t[gr_label_0 == 3, ]), cov(data_set_t))
    x <- rbind(x_1, x_2, x_3)
  }
  # matplot(t(x), type = "l", lty = 1, col = rep(1:3, n_vec))
  
  obs1 <- simul.obs(n_vec[1], grd = ti, d = 1.3 * max(ti), f = .5 * max(ti))
  obs2 <- simul.obs(n_vec[2], grd = ti, d = 1.3 * max(ti), f = .5 * max(ti))
  obs3 <- simul.obs(n_vec[3], grd = ti, d = 1.3 * max(ti), f = .5 * max(ti))
  # partially observed curves
  x1 <- ifelse(obs1, x[gr_label == 1, ], NA)
  x2 <- ifelse(obs2, x[gr_label == 2, ], NA)
  x3 <- ifelse(obs3, x[gr_label == 3, ], NA)
  
  # ANOVA for partially observed functional data
  # convert the data x1, x2, x3 into a data.table in long (tall) format
  x.wide <- data.table(id = 1:(sum(n_vec)), rbind(x1, x2, x3), 
                       grp = c(rep(1L, n_vec[1]), rep(2L, n_vec[2]), rep(3L, n_vec[3])))
  x <- melt(x.wide, id.vars = c("id", "grp"), value.name = "x")
  x[, variable := NULL]
  x[, ti := ti, by = "id,grp"]
  setcolorder(x, c("id", "grp", "ti", "x"))
  test_res <- meantest(x)
  
  if (is.null(test_res$message)) {
    results[i_nr, ] <- c(test_res$p.l2.test, test_res$p.proj.test, test_res$p.proj.test.boot)
  } else if (test_res$message == "no usable bootstrap samples") {
    temp <- test_res$message
    while (temp == "no usable bootstrap samples") {
      if (hipo == "H0") {
        x <- rmvnorm(sum(n_vec), colMeans(data_set_t), cov(data_set_t))
      } else if (hipo == "H1") {
        x_1 <- rmvnorm(n_vec[1], colMeans(data_set_t[gr_label_0 == 1, ]), cov(data_set_t))
        x_2 <- rmvnorm(n_vec[2], colMeans(data_set_t[gr_label_0 == 2, ]), cov(data_set_t))
        x_3 <- rmvnorm(n_vec[3], colMeans(data_set_t[gr_label_0 == 3, ]), cov(data_set_t))
        x <- rbind(x_1, x_2, x_3)
      }
      # matplot(t(x), type = "l", lty = 1, col = rep(1:3, n_vec))
      
      obs1 <- simul.obs(n_vec[1], grd = ti, d = 1.3 * max(ti), f = .5 * max(ti))
      obs2 <- simul.obs(n_vec[2], grd = ti, d = 1.3 * max(ti), f = .5 * max(ti))
      obs3 <- simul.obs(n_vec[3], grd = ti, d = 1.3 * max(ti), f = .5 * max(ti))
      # partially observed curves
      x1 <- ifelse(obs1, x[gr_label == 1, ], NA)
      x2 <- ifelse(obs2, x[gr_label == 2, ], NA)
      x3 <- ifelse(obs3, x[gr_label == 3, ], NA)
      
      # ANOVA for partially observed functional data
      # convert the data x1, x2, x3 into a data.table in long (tall) format
      x.wide <- data.table(id = 1:(sum(n_vec)), rbind(x1, x2, x3), 
                           grp = c(rep(1L, n_vec[1]), rep(2L, n_vec[2]), rep(3L, n_vec[3])))
      x <- melt(x.wide, id.vars = c("id", "grp"), value.name = "x")
      x[, variable := NULL]
      x[, ti := ti, by = "id,grp"]
      setcolorder(x, c("id", "grp", "ti", "x"))
      test_res <- meantest(x)
      if (is.null(test_res$message)) {
        temp <- "dobrze"
      }
    }
    results[i_nr, ] <- c(test_res$p.l2.test, test_res$p.proj.test, test_res$p.proj.test.boot)
  }
}
(results_final <- 100 * colMeans(results <= 0.05))
(time_total <- Sys.time() - time_1)

save(results_final, time_total,
     file = paste(paste("z_sim_rde_4", hipo, "d", delta, sep = "_"), "RData", sep = "."))
