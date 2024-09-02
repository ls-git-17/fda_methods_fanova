# 4. Multivariate analysis of variance

# for FMANOVA test by Gorecki and Smaga (2017)
library(fdANOVA)
# for Canadian weather data
library(fda)

# temperature and precipitation data
data_set_t <- t(CanadianWeather$dailyAv[,, "Temperature.C"])
data_set_p <- t(CanadianWeather$dailyAv[,, "Precipitation.mm"])
#  data for both features
data_set_t_p <- vector("list", 2)
data_set_t_p[[1]] <- t(data_set_t)
data_set_t_p[[2]] <- t(data_set_p)

# vector of group labels
gr_label <- rep(c(1, 2, 3), c(15, 15, 5))

# the tests based on a basis function representation
fmanova_ptbfr <- fmanova.ptbfr(data_set_t_p, gr_label, maxK = 31)
summary(fmanova_ptbfr)

# the tests based on random projections with the Gaussian white noise
set.seed(123)
fmanova_trp_1 <- fmanova.trp(data_set_t_p, gr_label, k = c(1, 5, 10, 15, 20))
summary(fmanova_trp_1)
plot(x = fmanova_trp_1)

# the tests based on random projections with the Brownian motion
set.seed(123)
fmanova_trp_2 <- fmanova.trp(data_set_t_p, gr_label, k = c(1, 5, 10, 15, 20), 
                             projection = "BM")
summary(fmanova_trp_2)
plot(x = fmanova_trp_2)

# implementation of the tests by Qiu, Chen, and Zhang (2021)

# stack or vec function
vec_fun <- function(x) {
  x_temp <- matrix(x[, 1], nrow(x))
  for (i in 2:ncol(x)) {
    x_temp <- rbind(x_temp, matrix(x[, i], nrow(x)))
  }
  return(x_temp)
}

# group mean and covariance calculation
# x - list of k elements corresponding to groups. Each element representing 
#     a group is a list of p elements corresponding to functional variables, 
#     and each such element representing a functional variable is a matrix
#     n_i times ntp (number of design time points) of descrete observations 
#     in design time points.
mean_cov_m <- function(x) {
  kk <- length(x)
  pp <- length(x[[1]])
  for (i in 2:kk) {
    if (pp != length(x[[i]])) stop("check the dimensions of functional data in 
                                   different groups")
  }
  ntp <- ncol(x[[1]][[1]])
  for (i in seq_len(kk)) {
    for (j in seq_len(pp)) {
      if (ntp != ncol(x[[i]][[j]])) stop("check the number of design time points 
                                         of functional data in different groups 
                                         and dimensions")
    }
  }
  n_vec <- numeric(kk)
  for (i in seq_len(kk)) n_vec[i] <- nrow(x[[i]][[1]])
  for (i in seq_len(kk)) {
    for (j in seq_len(pp)) {
      if (n_vec[i] != nrow(x[[i]][[j]])) stop("check the number of observations 
                                              of functional data in different groups")
    }
  }
  nn <- sum(n_vec)
  
  # sample vectors of mean functions
  gr_means <- vector("list", kk)
  for (i in seq_len(kk)) {
    mean_temp <- matrix(0, pp, ntp)
    for (j in 1:pp) {
      mean_temp[j, ] <- colMeans(x[[i]][[j]])
    }
    gr_means[[i]] <- mean_temp
  }
  gr_means_stack <- gr_means[[1]]
  for (i in 2:kk) {
    gr_means_stack <- rbind(gr_means_stack, gr_means[[i]])
  }
  
  x_obs <- vector("list", kk)
  for (i in seq_len(kk)) {
    x_temp <- vector("list", n_vec[i])
    for (j in seq_len(n_vec[i])) {
      x_temp[[j]] <- matrix(0, pp, ntp)
      for (l in seq_len(pp)) {
        x_temp[[j]][l, ] <- x[[i]][[l]][j, ] - gr_means[[i]][l, ]
      }
    }
    x_obs[[i]] <- x_temp
  }
  
  gr_cov <- vector("list", kk)
  for (i in seq_len(kk)) {
    cov_temp <- matrix(0, pp * ntp, pp * ntp)
    for (j in seq_len(n_vec[i])) {
      x_ij_stack <- vec_fun(x_obs[[i]][[j]])
      cov_temp <- cov_temp + x_ij_stack %*% t(x_ij_stack)
    }
    gr_cov[[i]] <- cov_temp / (n_vec[i] - 1)
  }
  
  return(list(nn = nn, ntp = ntp, kk = kk, n_vec = n_vec, pp = pp,
              gr_means = gr_means_stack, gr_cov = gr_cov, x_obs = x_obs))
}

# test statistics
qiu_chen_zhang_ts <- function(x) {
  mc <- mean_cov_m(x)
  gr_means <- mc$gr_means
  n_vec <- mc$n_vec
  gr_cov <- mc$gr_cov
  nn <- mc$nn
  ntp <- mc$ntp
  pp <- mc$pp
  x_obs <- mc$x_obs
  
  y_1_m_y_2 <- gr_means[1:pp, ] - gr_means[(pp + 1):(2 * pp), ]
  y_1_m_y_2_t <- t(y_1_m_y_2)
  
  gamma_hat <- matrix(0, nrow = nrow(gr_cov[[1]]), ncol = ncol(gr_cov[[1]]))
  for (i in seq_len(mc$kk)) {
    gamma_hat <- gamma_hat + (n_vec[i] - 1) * gr_cov[[i]]
  }
  gamma_hat <- gamma_hat / (sum(n_vec) - 2)
  
  t_point <- numeric(ntp)
  for (i in seq_len(ntp)) {
    t_point[i] <- (prod(n_vec) / nn) * y_1_m_y_2_t[i, ] %*% 
      solve(gamma_hat[(i * pp + 1 - pp):(i * pp), (i * pp + 1 - pp):(i * pp)]) %*% y_1_m_y_2[, i]
  }
  
  return(list(mean(t_point), max(t_point), nn, ntp, pp, n_vec, gr_means, gamma_hat, x_obs))
}

# tests
qiu_chen_zhang_tests <- function(x, a = 0, b = 1, n_boot = 1000) {
  temp <- qiu_chen_zhang_ts(x)
  nn <- temp[[3]]
  ntp <- temp[[4]]
  pp <- temp[[5]]
  n_vec <- temp[[6]]
  gr_means <- temp[[7]]
  gamma_hat <- temp[[8]]
  x_obs <- temp[[9]]
  vv <- c(x_obs[[1]], x_obs[[2]])
  
  # test based on the W-S approximation for integral statistic
  trace_st <- 0
  for (s in seq_len(ntp)) {
    for (t in seq_len(ntp)) {
      gamma_s <- gamma_hat[(s * pp + 1 - pp):(s * pp), (s * pp + 1 - pp):(s * pp)]
      gamma_t <- gamma_hat[(t * pp + 1 - pp):(t * pp), (t * pp + 1 - pp):(t * pp)]
      gamma_s_2 <- solve(eigen(gamma_s)$vectors %*% diag(sqrt(abs(eigen(gamma_s)$values))) %*% t(eigen(gamma_s)$vectors))
      gamma_t_2 <- solve(eigen(gamma_t)$vectors %*% diag(sqrt(abs(eigen(gamma_t)$values))) %*% t(eigen(gamma_t)$vectors))
      gamma_hat_s <- gamma_s_2 %*% gamma_hat[(s * pp + 1 - pp):(s * pp), (t * pp + 1 - pp):(t * pp)] %*% gamma_t_2
      
      trace_st <- trace_st + sum(diag(gamma_hat_s %*% t(gamma_hat_s)))
    }
  }
  trace_st <- trace_st / (ntp^2)
  beta_hat <- 1 / (pp * (b - a)) * trace_st
  d_hat <- (pp^2) * (b - a)^2 / trace_st
  p_val_1 <- 1 - pchisq(qiu_chen_zhang_ts(x)[[1]] / beta_hat, d_hat)
  
  # test based on bootstrap for supremum statistic
  sup_boot <- numeric(n_boot)
  for (i_b in seq_len(n_boot)) {
    vv_b <- vv[sample(1:nn, replace = TRUE)]
    vv_b_l <- list(vv_b[1:n_vec[1]], vv_b[-(1:n_vec[1])])
    vv_b_l_2 <- vector("list", 2)
    for (i_k in seq_len(2)) {
      vv_b_l_temp <- vector("list", pp)
      for (i_p in seq_len(pp)) {
        vv_b_l_temp[[i_p]] <- matrix(0, n_vec[i_k], ntp)
        for (i_n in seq_len(n_vec[i_k])) {
          vv_b_l_temp[[i_p]][i_n, ] <- vv_b_l[[i_k]][[i_n]][i_p, ]
        }
      }
      vv_b_l_2[[i_k]] <- vv_b_l_temp
    }
    sup_boot[i_b] <- qiu_chen_zhang_ts(vv_b_l_2)[[2]]
  }
  p_val_2 <- mean(sup_boot > temp[[2]])
  
  return(c(p_val_1, p_val_2))
}

# faster version using a part of code in C++

library(Rcpp)
sourceCpp("sec_4_af_gcz.cpp")

mean_cov_m_qcz_cpp <- function(x) {
  kk <- length(x)
  pp <- length(x[[1]])
  for (i in 2:kk) {
    if (pp != length(x[[i]])) stop("check the dimensions of functional data in different groups")
  }
  ntp <- ncol(x[[1]][[1]])
  for (i in seq_len(kk)) {
    for (j in seq_len(pp)) {
      if (ntp != ncol(x[[i]][[j]])) stop("check the number of design time points of functional data in different groups and dimensions")
    }
  }
  n_vec <- numeric(kk)
  for (i in seq_len(kk)) n_vec[i] <- nrow(x[[i]][[1]])
  for (i in seq_len(kk)) {
    for (j in seq_len(pp)) {
      if (n_vec[i] != nrow(x[[i]][[j]])) stop("check the number of observations of functional data in different groups")
    }
  }
  nn <- sum(n_vec)
  # sample vectors of mean functions
  gr_means <- vector("list", kk)
  for (i in seq_len(kk)) {
    mean_temp <- matrix(0, pp, ntp)
    for (j in 1:pp) {
      mean_temp[j, ] <- colMeans(x[[i]][[j]])
    }
    gr_means[[i]] <- mean_temp
  }
  gr_means_stack <- gr_means[[1]]
  for (i in 2:kk) {
    gr_means_stack <- rbind(gr_means_stack, gr_means[[i]])
  }
  
  temp_cpp_1 <- obs_cent_cpp(x, gr_means, kk, pp, ntp, n_vec)
  temp_cpp_2 <- gr_cov_cpp_2(temp_cpp_1, kk, pp, ntp, n_vec)
  
  return(list(nn = nn, ntp = ntp, kk = kk, n_vec = n_vec, pp = pp,
              gr_means = gr_means_stack, gr_cov = temp_cpp_2, x_obs = temp_cpp_1))
}

qiu_chen_zhang_ts_cpp <- function(x) {
  mc <- mean_cov_m_qcz_cpp(x)
  gr_means <- mc$gr_means
  n_vec <- mc$n_vec
  gr_cov <- mc$gr_cov
  nn <- mc$nn
  ntp <- mc$ntp
  pp <- mc$pp
  x_obs <- mc$x_obs
  
  y_1_m_y_2 <- gr_means[1:pp, ] - gr_means[(pp + 1):(2 * pp), ]
  y_1_m_y_2_t <- t(y_1_m_y_2)
  
  gamma_hat <- matrix(0, nrow = nrow(gr_cov[[1]]), ncol = ncol(gr_cov[[1]]))
  for (i in seq_len(mc$kk)) {
    gamma_hat <- gamma_hat + (n_vec[i] - 1) * gr_cov[[i]]
  }
  gamma_hat <- gamma_hat / (sum(n_vec) - 2)
  
  t_point <- numeric(ntp)
  for (i in seq_len(ntp)) {
    t_point[i] <- (prod(n_vec) / nn) * y_1_m_y_2_t[i, ] %*% 
      solve(gamma_hat[(i * pp + 1 - pp):(i * pp), (i * pp + 1 - pp):(i * pp)]) %*% y_1_m_y_2[, i]
  }
  
  return(list(mean(t_point), max(t_point), nn, ntp, pp, n_vec, gr_means, gamma_hat, x_obs))
}

qiu_chen_zhang_tests_cpp <- function(x, a = 0, b = 1, n_boot = 1000) {
  temp <- qiu_chen_zhang_ts_cpp(x)
  nn <- temp[[3]]
  ntp <- temp[[4]]
  pp <- temp[[5]]
  n_vec <- temp[[6]]
  gr_means <- temp[[7]]
  gamma_hat <- temp[[8]]
  x_obs <- temp[[9]]
  vv <- c(x_obs[[1]], x_obs[[2]])
  
  # test based on the W-S approximation for integral statistic
  trace_st <- 0
  for (s in seq_len(ntp)) {
    for (t in seq_len(ntp)) {
      gamma_s <- gamma_hat[(s * pp + 1 - pp):(s * pp), (s * pp + 1 - pp):(s * pp)]
      gamma_t <- gamma_hat[(t * pp + 1 - pp):(t * pp), (t * pp + 1 - pp):(t * pp)]
      gamma_s_2 <- solve(eigen(gamma_s)$vectors %*% diag(sqrt(abs(eigen(gamma_s)$values))) %*% t(eigen(gamma_s)$vectors))
      gamma_t_2 <- solve(eigen(gamma_t)$vectors %*% diag(sqrt(abs(eigen(gamma_t)$values))) %*% t(eigen(gamma_t)$vectors))
      gamma_hat_s <- gamma_s_2 %*% gamma_hat[(s * pp + 1 - pp):(s * pp), (t * pp + 1 - pp):(t * pp)] %*% gamma_t_2
      
      trace_st <- trace_st + sum(diag(gamma_hat_s %*% t(gamma_hat_s)))
    }
  }
  trace_st <- trace_st / (ntp^2)
  beta_hat <- 1 / (pp * (b - a)) * trace_st
  d_hat <- (pp^2) * (b - a)^2 / trace_st
  p_val_1 <- 1 - pchisq(qiu_chen_zhang_ts_cpp(x)[[1]] / beta_hat, d_hat)
  
  # test based on bootstrap for supremum statistic
  sup_boot <- numeric(n_boot)
  for (i_b in seq_len(n_boot)) {
    vv_b <- vv[sample(1:nn, replace = TRUE)]
    vv_b_l <- list(vv_b[1:n_vec[1]], vv_b[-(1:n_vec[1])])
    vv_b_l_2 <- vector("list", 2)
    for (i_k in seq_len(2)) {
      vv_b_l_temp <- vector("list", pp)
      for (i_p in seq_len(pp)) {
        vv_b_l_temp[[i_p]] <- matrix(0, n_vec[i_k], ntp)
        for (i_n in seq_len(n_vec[i_k])) {
          vv_b_l_temp[[i_p]][i_n, ] <- vv_b_l[[i_k]][[i_n]][i_p, ]
        }
      }
      vv_b_l_2[[i_k]] <- vv_b_l_temp
    }
    sup_boot[i_b] <- qiu_chen_zhang_ts_cpp(vv_b_l_2)[[2]]
  }
  p_val_2 <- mean(sup_boot > temp[[2]])
  
  return(c(p_val_1, p_val_2))
}

# application of the tests by Qiu, Chen, and Zhang (2021) for comparing 
# the temperature and precipitation in Eastern and Western Canada
data_set <- list(list(data_set_t[gr_label == 1, ], data_set_p[gr_label == 1, ]),
                 list(data_set_t[gr_label == 2, ], data_set_p[gr_label == 2, ]))
library(tictoc)
tic()
qiu_chen_zhang_tests(data_set)
toc()
## 0.0009174405 0.0000000000
## 401.401 sec elapsed
tic()
qiu_chen_zhang_tests_cpp(data_set)
toc()
## 0.0009174405 0.0000000000
## 215.814 sec elapsed
