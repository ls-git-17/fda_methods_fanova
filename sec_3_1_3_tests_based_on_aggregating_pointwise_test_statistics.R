# 3. Univariate analysis of variance
# 3.1.3. Tests based on aggregating pointwise test statistics

# for implementation of test
library(fdANOVA)
# for Canadian weather data
library(fda)

# temperature data
data_set_t <- t(CanadianWeather$dailyAv[,, "Temperature.C"])

# vector of group labels
gr_label <- rep(c(1, 2, 3), c(15, 15, 5))

# Table 2. Results of tests based on aggregating pointwise test statistics for 
# comparing temperature mean functions in three regions in Canada.

fanova_apts <- fanova.tests(t(data_set_t), gr_label, 
                            test = c("L2N", "L2B", "L2b", 
                                     "FN", "FB", "Fb", 
                                     "GPF", "Fmaxb"))
summary(fanova_apts)

# Figure 5. Pointwise test statistics

# function for pointwise test statistics
ssr_f_point <- function(x, group_label) {
  group_label_0 <- unique(group_label)
  l <- length(group_label_0)
  n_i <- numeric(l)
  for (i in 1:l) n_i[i] <- sum(group_label == group_label_0[i])
  n <- ncol(x)
  p <- nrow(x)
  if (n != length(group_label)) {
    stop("number of observations (number of columns in x) and number of elements in vector of group labels (group_label) must be the same")
  }
  mu0 <- rowMeans(x)
  SSR <- 0
  SSE <- 0
  for (i in seq_len(l)) {
    xi <- x[, group_label == group_label_0[i]]
    mui <- rowMeans(xi)
    zi <- t(xi) - as.matrix(rep(1, n_i[i])) %*% mui
    SSR <- SSR + n_i[i] * (mui - mu0)^2
    SSE <- SSE + colSums(zi^2)
  }
  return(list(SSR = SSR, F = (SSR / (l - 1)) / (SSE / (n - l))))
}
ssr_f <- ssr_f_point(t(data_set_t), gr_label)

# plot
par(mfrow = c(1, 2), mar = c(4, 2, 2, 0.1))
plot(ssr_f$SSR, type = "l", xlab = "Day", ylab = "", main = "SSR(t)")
plot(ssr_f$F, type = "l", xlab = "Day", ylab = "", main = "F(t)")
