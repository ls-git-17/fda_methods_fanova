# 5. Simulation studies based on real data examples
# 5.2 ANOVA for independent functional samples

hipo <- "H0" # H0, H1
delta <- 1

# for implementation of test
library(fdANOVA)
# for Canadian weather data
library(fda)
# for random number generators
library(mvtnorm)
# for implementation of graphical functional ANOVA test
library(GET)

# temperature data
data_set_t <- t(CanadianWeather$dailyAv[,, "Temperature.C"])
# str(data_set_t)

# vector of sample sizes
n_vec <- delta * c(15, 15, 5)

# vector of group labels
gr_label_0 <- rep(c(1, 2, 3), c(15, 15, 5))
gr_label <- rep(c(1, 2, 3), n_vec)

nr <- 1000
results <- matrix(NA, nr, 28)
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
  
  # FANOVA test based on basis expansion
  fanova <- fanova.tests(x = t(x), group.label = gr_label, test = "FP",
                         params = list(paramFP = list(maxK = 21)))
  # the tests based on random projections with the Gaussian white noise
  fanova_trp_1 <- fanova.tests(t(x), gr_label, test = "TRP",
                               params = list(paramTRP = list(k = c(10, 20, 30), B.TRP = 1000)))
  # the tests based on random projections with the Brownian motion
  fanova_trp_2 <- fanova.tests(t(x), gr_label, test = "TRP",
                               params = list(paramTRP = list(k = c(10, 20, 30), B.TRP = 1000, 
                                                             projection = "BM")))
  # tests based on aggregating pointwise test statistics
  fanova_apts <- fanova.tests(t(x), gr_label, 
                              test = c("L2N", "L2B", "L2b", 
                                       "FN", "FB", "Fb", 
                                       "GPF", "Fmaxb"))
  # graphical functional ANOVA test
  # vector of group labels
  gr_label_f_2 <- factor(rep(c("Eastern", "Western", "Northern"), n_vec),
                         levels = c("Eastern", "Western", "Northern"))
  cset <- create_curve_set(list(r = 1:365, obs = t(x)))
  res_gfanova <- graph.fanova(nsim = 1000, curve_set = cset,
                              groups = gr_label_f_2, variances = "equal",
                              contrasts = FALSE)
  
  results[i_nr, ] <- c(fanova$FP$pvalueFP,
                       fanova_trp_1$TRP$pvalues.anova, fanova_trp_1$TRP$pvalues.ATS, fanova_trp_1$TRP$pvalues.WTPS,
                       fanova_trp_2$TRP$pvalues.anova, fanova_trp_2$TRP$pvalues.ATS, fanova_trp_2$TRP$pvalues.WTPS,
                       fanova_apts$L2N$pvalueL2N, fanova_apts$L2B$pvalueL2B, fanova_apts$L2b$pvalueL2b,
                       fanova_apts$FN$pvalueFN, fanova_apts$FB$pvalueFB, fanova_apts$Fb$pvalueFb,
                       fanova_apts$GPF$pvalueGPF, fanova_apts$Fmaxb$pvalueFmaxb,
                       attr(res_gfanova, "p"))
}
(results_final <- 100 * colMeans(results <= 0.05))
(time_total <- Sys.time() - time_1)

save(results_final, time_total,
     file = paste(paste("z_sim_rde_3", hipo, "d", delta, sep = "_"), "RData", sep = "."))

# FANOVA test based on basis expansion
results_final[1]
# the tests based on random projections with the Gaussian white noise
results_final[2:10]
# the tests based on random projections with the Brownian motion
results_final[11:19]
# tests based on aggregating pointwise test statistics
# L2-norm-based tests
results_final[20:22]
# F-type tests
results_final[23:25]
# GPF and F_max tests
results_final[26:27]
# graphical functional ANOVA test
results_final[28]

plot(results_final, type = "l")
text(1:28, results_final, 1:28)
