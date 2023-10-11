# 3. Univariate analysis of variance
# 3.5. ANOVA for partially observed functional data

# for Canadian weather data
library(fda)
# for image.plot() function
library(fields)
# for data storing
library(data.table)

# download the supplementary code for the paper
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

# vector of group labels
gr_label <- rep(c(1, 2, 3), c(15, 15, 5))

# generate partially observed functional data from temperature data
ti <- 1:365
n_1 <- 15
n_2 <- 15
n_3 <- 5
set.seed(1234)
obs1 <- simul.obs(n_1, grd = ti, d = 1.3 * max(ti), f = .5 * max(ti))
obs2 <- simul.obs(n_2, grd = ti, d = 1.3 * max(ti), f = .5 * max(ti))
obs3 <- simul.obs(n_3, grd = ti, d = 1.3 * max(ti), f = .5 * max(ti))
# partially observed curves
x1 <- ifelse(obs1, data_set_t[gr_label == 1, ], NA)
x2 <- ifelse(obs2, data_set_t[gr_label == 2, ], NA)
x3 <- ifelse(obs3, data_set_t[gr_label == 3, ], NA)

# Figure 7. Partially observed temperature data and the summary of missingness
par(mfrow = c(1, 3), mar = c(4, 4, 2, 0.1))
matplot(ti, t(x2), type = "l", lty = 1, lwd = 0.5, 
        main = "Western Canada - partially observed curves", 
        xlab = "Day", ylab = "Temperature (C)", col = 1)
plot(ti, colMeans(!is.na(x2)), type = "l", 
     main = "Cross-sectional proportion\nof observed values", 
     xlab = "Day", ylab = "")
image.plot(ti, ti, crossprod(!is.na(x2)) / n_2, 
           main = "Proportion of pairwise\ncomplete observations",
           xlab = "Day", ylab = "Day")

# plot for all three groups of stations
par(mfrow = c(1, 3), mar = c(4, 4, 2, 0.1))
matplot(ti, t(x1), type = "l", lty = 1, main = "Eastern Canada", 
        xlab = "Day", ylab = "Temperature (C)", col = 1)
matplot(ti, t(x2), type = "l", lty = 1, main = "Western Canada", 
        xlab = "Day", ylab = "Temperature (C)", col = 1)
matplot(ti, t(x3), type = "l", lty = 1, main = "Northern Canada", 
        xlab = "Day", ylab = "Temperature (C)", col = 1)

# ANOVA for partially observed functional data

# convert the data x1,x2,x3 into a data.table in long (tall) format
x.wide <- data.table(id = 1:(n_1 + n_2 + n_3), rbind(x1, x2, x3), 
                     grp = c(rep(1L, n_1), rep(2L, n_2), rep(3L, n_3)))
x <- melt(x.wide, id.vars = c("id", "grp"), value.name = "x")
x[, variable := NULL]
x[, ti := ti, by = "id,grp"]
setcolorder(x, c("id", "grp", "ti", "x"))
set.seed(1)
meantest(x)
