# One-day VaR from ARMA(1,1)-GARCH(2,2), normal vs SGED, with Kupiec and
# Christoffersen backtests.
#
# Both models are fitted on the full sample and then backtested on it, so this
# is in-sample. A proper exercise would refit on a rolling window the way
# rolling_forecast.R does. The comparison between the two distributions still
# holds since both get the same advantage.

library(quantmod)
library(fGarch)

set.seed(123)

dir.create("output", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
sink("output/var_backtest.txt", split = TRUE)

invisible(getSymbols(Symbols = "^GSPC", from = "2000-01-01", to = "2022-11-12", auto.assign = TRUE))
log_returns <- na.omit(as.numeric(dailyReturn(GSPC, type = "log")))

norm_fit <- garchFit(~ arma(1, 1) + garch(2, 2), data = log_returns,
                     cond.dist = "norm", trace = FALSE)
sged_fit <- garchFit(~ arma(1, 1) + garch(2, 2), data = log_returns,
                     cond.dist = "sged", trace = FALSE)

shape <- sged_fit@fit$par[["shape"]]
skew <- sged_fit@fit$par[["skew"]]
cat("SGED shape:", round(shape, 4), " skew:", round(skew, 4), "\n\n")

# Conditional mean and sd from each fit
mu_norm <- as.numeric(fitted(norm_fit))
sd_norm <- as.numeric(volatility(norm_fit, type = "sigma"))
mu_sged <- as.numeric(fitted(sged_fit))
sd_sged <- as.numeric(volatility(sged_fit, type = "sigma"))

# Kupiec unconditional coverage, chi-square 1 df
kupiec <- function(hits, alpha) {
  n <- length(hits)
  x <- sum(hits)
  pi_hat <- x / n
  lr <- -2 * ((n - x) * log(1 - alpha) + x * log(alpha) -
              (n - x) * log(1 - pi_hat) - x * log(pi_hat))
  c(rate = pi_hat, stat = lr, p = 1 - pchisq(lr, df = 1))
}

# Christoffersen independence, chi-square 1 df
christoffersen <- function(hits) {
  prev <- hits[-length(hits)]
  curr <- hits[-1]
  n00 <- sum(prev == 0 & curr == 0)
  n01 <- sum(prev == 0 & curr == 1)
  n10 <- sum(prev == 1 & curr == 0)
  n11 <- sum(prev == 1 & curr == 1)

  pi01 <- n01 / (n00 + n01)
  pi11 <- n11 / (n10 + n11)
  pi_all <- (n01 + n11) / (n00 + n01 + n10 + n11)

  lr <- -2 * ((n00 + n10) * log(1 - pi_all) + (n01 + n11) * log(pi_all) -
              n00 * log(1 - pi01) - n01 * log(pi01) -
              n10 * log(1 - pi11) - n11 * log(pi11))
  c(stat = lr, p = 1 - pchisq(lr, df = 1))
}

show_result <- function(label, alpha, hits) {
  k <- kupiec(hits, alpha)
  ci <- christoffersen(hits)
  cc <- k[["stat"]] + ci[["stat"]]
  cat(label, "alpha =", alpha, "\n")
  cat("  exceedances:", sum(hits), "of", length(hits), "\n")
  cat("  rate:", round(k[["rate"]], 4), "expected:", alpha, "\n")
  cat("  Kupiec LR:", round(k[["stat"]], 3), "p:", round(k[["p"]], 4), "\n")
  cat("  Christoffersen LR:", round(ci[["stat"]], 3), "p:", round(ci[["p"]], 4), "\n")
  cat("  combined LR:", round(cc, 3), "p:", round(1 - pchisq(cc, df = 2), 4), "\n\n")
}

# 1% VaR
var_norm_01 <- mu_norm + sd_norm * qnorm(0.01)
var_sged_01 <- mu_sged + sd_sged * qsged(0.01, nu = shape, xi = skew)
show_result("normal", 0.01, as.integer(log_returns < var_norm_01))
show_result("sged  ", 0.01, as.integer(log_returns < var_sged_01))

# 5% VaR
var_norm_05 <- mu_norm + sd_norm * qnorm(0.05)
var_sged_05 <- mu_sged + sd_sged * qsged(0.05, nu = shape, xi = skew)
show_result("normal", 0.05, as.integer(log_returns < var_norm_05))
show_result("sged  ", 0.05, as.integer(log_returns < var_sged_05))

sink()
