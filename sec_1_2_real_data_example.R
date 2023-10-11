# 1.2. Real data example

# Canadian weather data (both features)

# for Canadian weather data
library(fda)
# for plots of functional data and sample mean functions
library(fdANOVA)
# for graphics
library(ggplot2)
# for multiple plots in one figure
library(gridExtra)
# map of Canada
library(rworldmap)

# temperature and precipitation data
data_set_t <- t(CanadianWeather$dailyAv[,, "Temperature.C"])
data_set_p <- t(CanadianWeather$dailyAv[,, "Precipitation.mm"])

# vectors of group labels
gr_label <- rep(c(1, 2, 3), c(15, 15, 5))
gr_label_f_2 <- factor(rep(c("Eastern", "Western", "Northern"), c(15, 15, 5)),
                       levels = c("Eastern", "Western", "Northern"))

# Figure 1. Map of Canada with noted weather stations in three different regions
colors <- 1:3
names(colors) <- c("Eastern", "Western", "Northern")
cs <- factor(colors[rep(c("Eastern", "Western", "Northern"), c(15, 15, 5))])
can_data <- with(CanadianWeather, 
                 data.frame(long = -coordinates[, 2], 
                            lat = coordinates[, 1], 
                            group = 1, 
                            cs = cs))
canada_map <- map_data(map = "world", region = "Canada")
ggplot(canada_map, aes(x = long, y = lat, group = group)) + 
  geom_polygon(color = "black", fill = "white") +
  geom_point(data = can_data, aes(color = cs, shape = cs), size = 4) +
  scale_color_manual("Canada region", values = 2:4, labels = names(colors)) +
  scale_shape_manual("Canada region", values = 15:17, labels = names(colors)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

# Figure 2. Trajectories of the temperature and precipitation noted in Canadian 
# weather stations and the sample mean functions
par(mfrow = c(1, 4), mar = c(4, 4, 2, 0.1))
matplot(t(data_set_t[gr_label == 1, ]),
        type = "l", col = 1, lty = 1, ylim = c(-33, 22),
        main = "Eastern Canada", xlab = "Day", ylab = "Temperature (C)")
matplot(t(data_set_t[gr_label == 2, ]),
        type = "l", col = 1, lty = 1, ylim = c(-33, 22),
        main = "Western Canada", xlab = "Day", ylab = "Temperature (C)")
matplot(t(data_set_t[gr_label == 3, ]),
        type = "l", col = 1, lty = 1, ylim = c(-33, 22),
        main = "Northern Canada", xlab = "Day", ylab = "Temperature (C)")
matplot(matrix(c(colMeans(data_set_t[gr_label == 1, ]),
                 colMeans(data_set_t[gr_label == 2, ]),
                 colMeans(data_set_t[gr_label == 3, ])),
               ncol = 3),
        type = "l", col = 1, lty = 1:3, ylim = c(-33, 22),
        main = "Sample mean functions", xlab = "Day", ylab = "Temperature (C)")
legend("topright", legend = c("Eastern", "Western", "Northern"), lty = 1:3)

par(mfrow = c(1, 4), mar = c(4, 4, 2, 0.1))
matplot(t(data_set_p[gr_label == 1, ]),
        type = "l", col = 1, lty = 1, ylim = c(0, 17),
        main = "Eastern Canada", xlab = "Day", ylab = "Precipitation (mm)")
matplot(t(data_set_p[gr_label == 2, ]),
        type = "l", col = 1, lty = 1, ylim = c(0, 17),
        main = "Western Canada", xlab = "Day", ylab = "Precipitation (mm)")
matplot(t(data_set_p[gr_label == 3, ]),
        type = "l", col = 1, lty = 1, ylim = c(0, 17),
        main = "Northern Canada", xlab = "Day", ylab = "Precipitation (mm)")
matplot(matrix(c(colMeans(data_set_p[gr_label == 1, ]),
                 colMeans(data_set_p[gr_label == 2, ]),
                 colMeans(data_set_p[gr_label == 3, ])),
               ncol = 3),
        type = "l", col = 1, lty = 1:3, ylim = c(0, 17),
        main = "Sample mean functions", xlab = "Day", ylab = "Precipitation (mm)")
legend("topright", legend = c("Eastern", "Western", "Northern"), lty = 1:3)

# The same as in Figure 2, but with using ggplot2 and fdANOVA packages
plot_1 <- plotFANOVA(x = t(data_set_t[gr_label == 1, ])) + 
  ggtitle("Eastern Canada") + xlab("Day") + ylab("Temperature (C)") + ylim(-40, 22)
plot_2 <- plotFANOVA(x = t(data_set_t[gr_label == 2, ])) + 
  ggtitle("Western Canada") + xlab("Day") + ylab("Temperature (C)") + ylim(-40, 22)
plot_3 <- plotFANOVA(x = t(data_set_t[gr_label == 3, ])) + 
  ggtitle("Northern Canada") + xlab("Day") + ylab("Temperature (C)") + ylim(-40, 22)
plot_4 <- plotFANOVA(x = t(data_set_t), group.label = gr_label_f_2,
                     means = TRUE) +
  scale_color_manual("Canada region", values = 2:4, labels = names(colors)) +
  scale_shape_manual("Canada region", values = 15:17, labels = names(colors)) +
  ggtitle("Sample mean functions") + xlab("Day") + ylab("Temperature (C)") + ylim(-40, 22)
grid.arrange(plot_1, plot_2, plot_3, plot_4, ncol = 4)

plot_1 <- plotFANOVA(x = t(data_set_p[gr_label == 1, ])) + 
  ggtitle("Eastern Canada") + xlab("Day") + ylab("Precipitation (mm)") + ylim(-0, 15)
plot_2 <- plotFANOVA(x = t(data_set_p[gr_label == 2, ])) + 
  ggtitle("Western Canada") + xlab("Day") + ylab("Precipitation (mm)") + ylim(-0, 15)
plot_3 <- plotFANOVA(x = t(data_set_p[gr_label == 3, ])) + 
  ggtitle("Northern Canada") + xlab("Day") + ylab("Precipitation (mm)") + ylim(-0, 15)
plot_4 <- plotFANOVA(x = t(data_set_p), group.label = gr_label_f_2,
                     means = TRUE) +
  scale_color_manual("Canada region", values = 2:4, labels = names(colors)) +
  scale_shape_manual("Canada region", values = 15:17, labels = names(colors)) +
  ggtitle("Sample mean functions") + xlab("Day") + ylab("Precipitation (mm)") + ylim(-0, 15)
grid.arrange(plot_1, plot_2, plot_3, plot_4, ncol = 4)
