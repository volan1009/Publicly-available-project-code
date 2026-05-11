## install.packages("nlme")
library(nlme)

set.seed(2026)

## ---- Parameters ----
wmft_baseline <- 33         # baseline WMFT (seconds)
log_wmft0     <- log(wmft_baseline)
beta_trt      <- -0.05      # log-scale effect per treated 2-week period
phi           <- 0.20       # AR(1) coefficient for level fluctuation
sigma_level2  <- 0.01       # innovation variance of the level fluctuation
sigma_obs2    <- 0          # measurement error variance
level_resid0  <- 0          # initial level residual

## ---- Design: 6 months with 2-week periods -> 12 periods ----
n_periods  <- 12
trt_period <- rep(c(0, 1), times = 6)    # ABAB... over 12 periods

weeks <- (0:n_periods) * 2

## ============================================================
## Single trajectory simulation (one subject)
## ============================================================

cum_trt <- c(0, cumsum(trt_period))

mean_log <- log_wmft0 + beta_trt * cum_trt

## AR(1) level fluctuation
level_innov <- rnorm(n_periods + 1, mean = 0, sd = sqrt(sigma_level2))
level_resid <- numeric(n_periods + 1)
level_resid[1] <- level_resid0
if (n_periods + 1 >= 2) {
  for (t in 2:(n_periods + 1)) {
    level_resid[t] <- phi * level_resid[t - 1] + level_innov[t]
  }
}

obs_noise <- rnorm(n_periods + 1, mean = 0, sd = sqrt(sigma_obs2))

wmft_det <- exp(mean_log)
wmft_obs <- exp(mean_log + level_resid + obs_noise)
Y_log <- mean_log + level_resid + obs_noise

df <- data.frame(
  period = 0:n_periods,
  weeks  = weeks,
  A      = c(NA, trt_period),
  CumTrt = cum_trt,
  mean_log = mean_log,
  level_resid = level_resid,
  wmft_det = wmft_det,
  Y_log   = Y_log,
  WMFT_s  = wmft_obs
)

## Plot single trajectory
plot(df$weeks, df$WMFT_s, type = "o", pch = 16, lty = 1, lwd = 1.5,
  xaxt = "n",
     xlab = "Weeks since baseline",
     ylab = "SWMFT-C task completion time (seconds)",
     main = "One time simulation trajectory of Willie Makeit")
axis(1, at = weeks, labels = weeks)
lines(df$weeks, df$wmft_det, type = "o", pch = 1, lty = 2, lwd = 1.5)
abline(v = df$weeks, lty = 3, col = "gray80")
text(df$weeks, df$WMFT_s,
  labels = sprintf("%.1f", df$WMFT_s), pos = 1, cex = 1)
text(df$weeks, df$wmft_det,
  labels = sprintf("%.1f", df$wmft_det), pos = 3, cex = 1)
legend("topright",
       legend = c("Observed trajectory", "Deterministic trajectory"),
       lty = c(1, 2), pch = c(16, 1), lwd = 1.5, bty = "n")

## GLS inference
fit_gls <- gls(
  Y_log ~ CumTrt,
  data = df,
  correlation = corAR1(form = ~ weeks)
)

sm <- summary(fit_gls)
print(sm)

beta_hat <- sm$tTable["CumTrt", "Value"]
se_hat   <- sm$tTable["CumTrt", "Std.Error"]

ci_low  <- beta_hat - 1.96 * se_hat
ci_high <- beta_hat + 1.96 * se_hat

c(beta_hat = beta_hat, ci_low = ci_low, ci_high = ci_high)

ratio_hat <- exp(beta_hat)
ratio_ci  <- exp(c(ci_low, ci_high))

c(ratio_hat = ratio_hat, ratio_ci_low = ratio_ci[1], ratio_ci_high = ratio_ci[2])


## ============================================================
## Monte Carlo: 1000 simulated trajectories + mean + 95% bands

n_sim <- 1000
wmft_mat <- matrix(NA, nrow = n_sim, ncol = n_periods + 1)

for (i in 1:n_sim) {
  level_innov <- rnorm(n_periods + 1, 0, sqrt(sigma_level2))
  level_resid <- numeric(n_periods + 1)
  level_resid[1] <- level_resid0
  for (t in 2:(n_periods + 1)) {
    level_resid[t] <- phi * level_resid[t - 1] + level_innov[t]
  }

  obs_noise <- rnorm(n_periods + 1, 0, sqrt(sigma_obs2))
  wmft_mat[i, ] <- exp(mean_log + level_resid + obs_noise)
}

mean_traj <- colMeans(wmft_mat)
lower_95 <- apply(wmft_mat, 2, quantile, 0.025)
upper_95 <- apply(wmft_mat, 2, quantile, 0.975)

plot(weeks, wmft_mat[1, ], type = "n",
     ylim = range(wmft_mat),
  xaxt = "n",
     xlab = "Weeks since baseline",
     ylab = "SWMFT-C task completion time (seconds)",
     main = "1000 time simulation trajectories of Willie Makeit")
axis(1, at = weeks, labels = weeks)

for (i in 1:n_sim) {
  lines(weeks, wmft_mat[i, ], col = rgb(0, 0, 0, 0.04), lwd = 0.7)
}

lines(weeks, exp(mean_log), col = "black", lwd = 1.5, lty = 2)
lines(weeks, mean_traj, col = "black", lwd = 2)
lines(weeks, lower_95, col = "black", lwd = 1.5, lty = 3)
lines(weeks, upper_95, col = "black", lwd = 1.5, lty = 3)
points(weeks, mean_traj, pch = 16, cex = 0.6)
points(weeks, lower_95, pch = 1, cex = 0.6)
points(weeks, upper_95, pch = 1, cex = 0.6)
text(weeks, mean_traj,
  labels = sprintf("%.1f", mean_traj), pos = 3, cex = 0.8)
text(weeks, lower_95,
  labels = sprintf("%.1f", lower_95), pos = 1, cex = 0.8)
text(weeks, upper_95,
  labels = sprintf("%.1f", upper_95), pos = 3, cex = 0.8)
abline(v = weeks, lty = 3, col = "gray80")
legend("topright",
       legend = c("Monte Carlo mean trajectory", "Deterministic trajectory", "95% empirical band"),
       lty = c(1, 2, 3), lwd = c(2, 1.5, 1.5), bty = "n")







