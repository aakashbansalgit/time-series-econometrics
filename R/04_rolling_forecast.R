# Rolling one-step variance forecasts, scored out of sample. AIC could not
# separate the top two specifications, so this checks whether either forecasts
# better. Takes a few minutes.

library(quantmod)
library(fGarch)

set.seed(123)

dir.create("output", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
sink("output/rolling_forecast.txt", split = TRUE)

invisible(getSymbols(Symbols = "^GSPC", from = "2000-01-01", to = "2022-11-12", auto.assign = TRUE))
log_returns <- na.omit(as.numeric(dailyReturn(GSPC, type = "log")))

n <- length(log_returns)
start_window <- 2000   # first forecast made after this many observations
refit_every <- 250     # refit annually, roll parameters between

specs <- list(
  "arma(1,1)+garch(2,2)" = ~ arma(1, 1) + garch(2, 2),
  "arma(0,2)+garch(2,2)" = ~ arma(0, 2) + garch(2, 2),
  "arma(0,0)+garch(1,1)" = ~ arma(0, 0) + garch(1, 1)
)

# Expanding window
roll_fcst <- function(formula) {
  pred <- rep(NA_real_, n)
  fit <- NULL
  for (t in start_window:(n - 1)) {
    if ((t - start_window) %% refit_every == 0 || is.null(fit)) {
      fit <- try(garchFit(formula, data = log_returns[1:t], trace = FALSE),
                 silent = TRUE)
      if (inherits(fit, "try-error")) {
        fit <- NULL
        next
      }
    }
    fc <- try(predict(fit, n.ahead = 1), silent = TRUE)
    if (!inherits(fc, "try-error")) {
      pred[t + 1] <- fc$standardDeviation[1]^2
    }
  }
  pred
}

# Squared return as the realised variance proxy
realised <- log_returns^2

qlike <- function(pred, actual) {
  ok <- !is.na(pred) & pred > 0
  mean(log(pred[ok]) + actual[ok] / pred[ok])
}

mse <- function(pred, actual) {
  ok <- !is.na(pred)
  mean((actual[ok] - pred[ok])^2)
}

results <- list()
for (nm in names(specs)) {
  cat("Running", nm, "\n")
  pred <- roll_fcst(specs[[nm]])
  results[[nm]] <- list(n = sum(!is.na(pred)),
                        mse = mse(pred, realised),
                        qlike = qlike(pred, realised))
}

cat("\nForecasts start at observation", start_window + 1, "of", n, "\n\n")
for (nm in names(results)) {
  r <- results[[nm]]
  cat(nm, "\n")
  cat("  forecasts:", r$n, "\n")
  cat("  MSE:", format(r$mse, digits = 4), "\n")
  cat("  QLIKE:", round(r$qlike, 5), "\n")
}

cat("\nLowest QLIKE:", names(which.min(sapply(results, function(r) r$qlike))), "\n")
cat("Lowest MSE:", names(which.min(sapply(results, function(r) r$mse))), "\n")

sink()
