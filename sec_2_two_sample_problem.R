# 2. Two sample problem

# for implementation of test
library(fda.usc)
# for Canadian weather data
library(fda)

# temperature data
data_set_t <- t(CanadianWeather$dailyAv[,, "Temperature.C"])

# vector of group labels
gr_label <- rep(c(1, 2, 3), c(15, 15, 5))

# data for Eastern, Western, and Northern Canada
ec <- data_set_t[gr_label == 1, ]
wc <- data_set_t[gr_label == 2, ]
nc <- data_set_t[gr_label == 3, ]

# Table 1. Results of two sample test for comparing temperature mean 
# functions in pairs of regions in Canada

# two sample test for Eastern-Western
res_1 <- fmean.test.fdata(ec, wc, method = "X2", npc = -0.95)
# test statistic
res_1$stat[1]
# p-value
res_1$pvalue[1]
# number of chosen functional principal components
res_1$p

# two sample test for Eastern-Northern
res_2 <- fmean.test.fdata(ec, nc, method = "X2", npc = -0.95)
# test statistic
res_2$stat[1]
# p-value
res_2$pvalue[1]
# number of chosen functional principal components
res_2$p

# two sample test for Western-Northern
res_3 <- fmean.test.fdata(wc, nc, method = "X2", npc = -0.95)
# test statistic
res_3$stat[1]
# p-value
res_3$pvalue[1]
# number of chosen functional principal components
res_3$p
