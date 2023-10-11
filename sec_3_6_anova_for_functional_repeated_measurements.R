# 3. Univariate analysis of variance
# 3.6. ANOVA for functional repeated measurements

# for implementation of test based on aggregating the pointwise test statistics
library(rmfanova)
# for DTI data set
library(refund)

# function for random projection tests by Smaga (2019b)
# x - list of matrices n times p corresponding to different objects (e.g., visits)
rpt_rfanova <- function(x, k = 30, projection = c("GAUSS", "BM")) {
  n <- nrow(x[[1]])
  p <- ncol(x[[1]])
  l <- length(x)
  p_values <- numeric(length(k))
  iik <- 0
  modulo <- function (z) { sqrt(sum(z^2)) }
  for (ik in k) {
    if (projection == "GAUSS") {
      z <- matrix(rnorm(p * ik), nrow = ik, ncol = p)
      modu <- apply(z, 1, modulo)
      z <- z / modu
      t_p <- numeric(ik)
      for (j in 1:ik) {
        temp <- c()
        for (jj in 1:l) {
          temp <- c(temp, c(x[[jj]] %*% z[j, ]))
        }
        data.temp <- data.frame(temp = temp, 
                                subject = as.factor(rep(1:n, l)), 
                                group = as.factor(rep(1:l, each = n)))
        model <- aov(temp ~ group + Error(subject / group), data = data.temp)
        t_p[j] <- summary(model)[[2]][[1]][[5]][1]
      }
      iik <- iik + 1
      p_values[iik] <- min(ik * t_p[order(t_p)] / 1:ik)
    } else {
      t_p <- numeric(ik)
      for (j in 1:ik) {
        bm.p <- cumsum(rnorm(p, mean = 0, sd = 1)) / sqrt(p)
        bm.p <- bm.p / modulo(bm.p)
        temp <- c()
        for (jj in 1:l) {
          temp <- c(temp, c(x[[jj]] %*% as.matrix(bm.p)))
        }
        data.temp <- data.frame(temp = temp, 
                                subject = as.factor(rep(1:n, l)), 
                                group = as.factor(rep(1:l, each = n)))
        model <- aov(temp ~ group + Error(subject / group), data = data.temp)
        t_p[j] <- summary(model)[[2]][[1]][[5]][1]
      }
      iik <- iik + 1
      p_values[iik] <- min(ik * t_p[order(t_p)] / 1:ik)
    }
  }
  return(list(p_values = p_values))
}

# preparation of the DTI data set, for details see Kurylo and Smaga (2023)
data(DTI)
# MS patients
DTI_ms <- DTI[DTI$case == 1, ]
miss_data <- c()
for (i in 1:340) if (any(is.na(DTI_ms$cca[i, ]))) miss_data <- c(miss_data, i)
DTI_ms <- DTI_ms[-miss_data, ]

DTI_ms_2 <- DTI_ms[DTI_ms$Nscans == 4, ]
xx <- vector("list", 4)
for (i in 1:4) {
  xx[[i]] <- DTI_ms_2$cca[DTI_ms_2$visit == i, ]
}
xx[[1]] <- xx[[1]][-14, ]
xx[[3]] <- xx[[3]][-14, ]
yy <- xx
for (i in seq_len(4)) yy[[i]] <- yy[[i]][1:17, ]

# Figure 8. Trajectories for the FA profiles at four visits
oldpar <- par(mfrow = c(1, 4), mar = c(4, 4, 4, 0.1))
matplot(t(yy[[1]]), type = "l", col = 1, lty = 1, xlab = "t", ylab = "FA",
        main = "Visit 1", xaxt = "n", ylim = c(0.29, 0.73))
axis(1, c(1, 15, 30, 45, 60, 75, 93), labels = c(1, 15, 30, 45, 60, 75, 93))
matplot(t(yy[[2]]), type = "l", col = 1, lty = 1, xlab = "t", ylab = "FA",
        main = "Visit 2", xaxt = "n", ylim = c(0.29, 0.73))
axis(1, c(1, 15, 30, 45, 60, 75, 93), labels = c(1, 15, 30, 45, 60, 75, 93))
matplot(t(yy[[3]]), type = "l", col = 1, lty = 1, xlab = "t", ylab = "FA",
        main = "Visit 3", xaxt = "n", ylim = c(0.29, 0.73))
axis(1, c(1, 15, 30, 45, 60, 75, 93), labels = c(1, 15, 30, 45, 60, 75, 93))
matplot(t(yy[[4]]), type = "l", col = 1, lty = 1, xlab = "t", ylab = "FA",
        main = "Visit 4", xaxt = "n", ylim = c(0.29, 0.73))
axis(1, c(1, 15, 30, 45, 60, 75, 93), labels = c(1, 15, 30, 45, 60, 75, 93))
par(oldpar)

# Figure 9. Sample mean functions of the FA at different visits and 
# pointwise test statistics
par(mfrow = c(1, 3), mar = c(4, 4, 2, 0.1))
# sample mean functions
pointwise_sample_mean_fun(yy, values = FALSE,
                          col = 1, lty = 1:4, xlab = "t", ylab = "FA", xaxt = "n")
axis(1, c(1, 15, 30, 45, 60, 75, 93), labels = c(1, 15, 30, 45, 60, 75, 93))
legend(x = 36, y = 0.64, legend = 1:4, col = 1, lty = 1:4, title = "Visit")
# pointwise SSA and F-type test statistics
par(mar = c(4, 2, 2, 0.1))
pointwise_ssa_test_statistic(yy, xlab = "t", xaxt = "n")
axis(1, c(1, 15, 30, 45, 60, 75, 93), labels = c(1, 15, 30, 45, 60, 75, 93))
pointwise_f_test_statistic(yy, xlab = "t", xaxt = "n")
axis(1, c(1, 15, 30, 45, 60, 75, 93), labels = c(1, 15, 30, 45, 60, 75, 93))

# Table 3. Results of the aggregating the pointwise test statistics procedures 
# by Kurylo and Smaga (2023) for the DTI data set
res <- rmfanova(yy)
summary(res, digits = 3)

# Figure 10. P-values of the random projection tests by Smaga (2019b) 
# for the DTI data set
set.seed(123)
results_gauss <- rpt_rfanova(yy, k = 5:30, projection = "GAUSS")
set.seed(123)
results_bm <- rpt_rfanova(yy, k = 5:30, projection = "BM")
par(mfrow = c(1, 2), mar = c(4, 4, 2, 0.1))
plot(5:30, results_gauss$p_values, type = "l", 
     main = "Gaussian white noise", xlab = "k", ylab = "p-value")
plot(5:30, results_bm$p_values, type = "l", 
     main = "Brownian motion", xlab = "k", ylab = "p-value")
