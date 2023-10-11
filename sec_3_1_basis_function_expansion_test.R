# 3. Univariate analysis of variance
# 3.1. Basis function expansion test

# for implementation of test
library(fdANOVA)
# for Canadian weather data
library(fda)

# temperature data
data_set_t <- t(CanadianWeather$dailyAv[,, "Temperature.C"])

# vector of group labels
gr_label <- rep(c(1, 2, 3), c(15, 15, 5))

# FANOVA test based on basis expansion
fanova <- fanova.tests(x = t(data_set_t), group.label = gr_label, test = "FP",
                       params = list(paramFP = list(maxK = 21)))
fanova
summary(fanova)

# Figure 3. In the left panel, the raw Canada temperature data are presented, 
# while in the right panel, we have the smoothed data in the Fourier basis 
# with seven basis functions.

# Fourier basis and its representation of data
fbasis <- create.fourier.basis(rangeval = c(0, 365), 7)
basis_exp <- Data2fd(1:365, t(data_set_t), fbasis)

# plot
par(mfrow = c(1, 2), mar = c(4, 4, 2, 0.1))
matplot(t(data_set_t), type = "l", col = 1, lty = gr_label,
        xlab = "Day", ylab = "Temperature (C)",
        main = "Canadian weather data set - raw data")
legend("bottom", legend = c("Eastern Canada", "Western Canada", "Northern Canada"),
       col = 1, lty = 1:3)
plot(basis_exp, col = 1, lty = gr_label, 
     xlab = "Day", ylab = "Temperature (C)",
     main = "Canadian weather data set - smoothed data")
legend("bottom", legend = c("Eastern Canada", "Western Canada", "Northern Canada"),
       col = 1, lty = 1:3)

# the same as above but with different colors for groups
par(mfrow = c(1, 2), mar = c(4, 4, 2, 0.1))
matplot(t(data_set_t), type = "l", col = gr_label, lty = 1,
        xlab = "Day", ylab = "Temperature (C)",
        main = "Canadian weather data set - raw data")
legend("bottom", legend = c("Eastern Canada", "Western Canada", "Northern Canada"),
       col = 1:3, lty = 1)
plot(basis_exp, col = gr_label, lty = 1, 
     xlab = "Day", ylab = "Temperature (C)",
     main = "Canadian weather data set - smoothed data")
legend("bottom", legend = c("Eastern Canada", "Western Canada", "Northern Canada"),
       col = 1:3, lty = 1)

# estimated coefficients of basis expansion of the data set
basis_exp$coefs
