# 3. Univariate analysis of variance
# 3.2. Tests based on random projections

# for implementation of test
library(fdANOVA)
# for Canadian weather data
library(fda)

# temperature data
data_set_t <- t(CanadianWeather$dailyAv[,, "Temperature.C"])

# vector of group labels
gr_label <- rep(c(1, 2, 3), c(15, 15, 5))

# the tests based on random projections with the Gaussian white noise
set.seed(123)
fanova_trp_1 <- fanova.tests(t(data_set_t), gr_label, test = "TRP",
                             params = list(paramTRP = list(k = 5:30, B.TRP = 1000)))
summary(fanova_trp_1)

# the tests based on random projections with the Brownian motion
set.seed(123)
fanova_trp_2 <- fanova.tests(t(data_set_t), gr_label, test = "TRP",
                             params = list(paramTRP = list(k = 5:30, B.TRP = 1000, 
                                                           projection = "BM")))
summary(fanova_trp_2)

# Figure 4. P-values of random projection test for comparing temperature 
# mean functions in three regions in Canada. k denotes the number 
# of projections used.

par(mfrow = c(1, 2), mar = c(4, 4, 2, 0.1))
matplot(x = 5:30, y = matrix(c(fanova_trp_1$TRP$pvalues.anova,
                               fanova_trp_1$TRP$pvalues.ATS,
                               fanova_trp_1$TRP$pvalues.WTPS),
                             ncol = 3),
        type = "l", col = 1, lty = 1:3,
        main = "Gaussian white noise", xlab = "k", ylab = "p-value")
legend("topright", legend = c("ANOVA", "ATS", "WTPS"), lty = 1:3)
matplot(x = 5:30, y = matrix(c(fanova_trp_2$TRP$pvalues.anova,
                               fanova_trp_2$TRP$pvalues.ATS,
                               fanova_trp_2$TRP$pvalues.WTPS),
                             ncol = 3),
        type = "l", col = 1, lty = 1:3,
        main = "Brownian motion", xlab = "k", ylab = "p-value")
legend("topright", legend = c("ANOVA", "ATS", "WTPS"), lty = 1:3)

# the same as above but using appropriate function from fdANOVA package
plot(fanova_trp_1)
plot(fanova_trp_2)
