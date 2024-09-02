# 3. Univariate analysis of variance
# 3.1.4. Graphical functional ANOVA test

# for implementation of test
library(GET)
# for Canadian weather data
library(fda)
# for graphics
library(ggplot2)

# temperature data
data_set_t <- t(CanadianWeather$dailyAv[,, "Temperature.C"])

# vector of group labels
gr_label_f_2 <- factor(rep(c("Eastern", "Western", "Northern"), c(15, 15, 5)),
                       levels = c("Eastern", "Western", "Northern"))

# graphical functional ANOVA test
cset <- create_curve_set(list(r = 1:365, obs = t(data_set_t)))
set.seed(123)
res_gfanova <- graph.fanova(nsim = 1000, curve_set = cset,
                            groups = gr_label_f_2, variances = "equal",
                            contrasts = FALSE)
res_gfanova

# Figure 6. The one-way graphical functional ANOVA test for equality of means 
# of the temperature in Canada in the three groups using the group means
plot(res_gfanova) + ggplot2::labs(x = "Day", y = "Temperature (C)")

# the same as above but black-white version for paper
par(mfrow = c(1, 3), mar = c(4, 4, 2, 0.1))

matplot(matrix(c(res_gfanova$Eastern$obs,
                 res_gfanova$Eastern$central,
                 res_gfanova$Eastern$lo,
                 res_gfanova$Eastern$hi),
               ncol = 4),
        type = "n", main = "Eastern", xlab = "Day", ylab = "Temperature (C)")
polygon(x = c(1:365, 365:1),
        y = c(res_gfanova$Eastern$lo, rev(res_gfanova$Eastern$hi)),
        border = NA, col = "grey")
lines(1:365, res_gfanova$Eastern$obs)
points((1:365)[res_gfanova$Eastern$obs > res_gfanova$Eastern$hi], 
       res_gfanova$Eastern$obs[res_gfanova$Eastern$obs > res_gfanova$Eastern$hi], 
       pch = 16)
lines(1:365, res_gfanova$Eastern$central, lty = 2)

matplot(matrix(c(res_gfanova$Western$obs,
                 res_gfanova$Western$central,
                 res_gfanova$Western$lo,
                 res_gfanova$Western$hi),
               ncol = 4),
        type = "n", main = "Western", xlab = "Day", ylab = "Temperature (C)")
polygon(x = c(1:365, 365:1),
        y = c(res_gfanova$Western$lo, rev(res_gfanova$Western$hi)),
        border = NA, col = "grey")
lines(1:365, res_gfanova$Western$obs)
points((1:365)[res_gfanova$Western$obs > res_gfanova$Western$hi], 
       res_gfanova$Western$obs[res_gfanova$Western$obs > res_gfanova$Western$hi], 
       pch = 16)
lines(1:365, res_gfanova$Western$central, lty = 2)

matplot(matrix(c(res_gfanova$Northern$obs,
                 res_gfanova$Northern$central,
                 res_gfanova$Northern$lo,
                 res_gfanova$Northern$hi),
               ncol = 4),
        type = "n", main = "Northern", xlab = "Day", ylab = "Temperature (C)")
polygon(x = c(1:365, 365:1),
        y = c(res_gfanova$Northern$lo, rev(res_gfanova$Northern$hi)),
        border = NA, col = "grey")
lines(1:365, res_gfanova$Northern$obs)
points((1:365)[res_gfanova$Northern$obs < res_gfanova$Northern$lo], 
       res_gfanova$Northern$obs[res_gfanova$Northern$obs < res_gfanova$Northern$lo], 
       pch = 16)
lines(1:365, res_gfanova$Northern$central, lty = 2)
