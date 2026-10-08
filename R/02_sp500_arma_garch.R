# Fits ARMA-GARCH models to S&P 500 daily log returns and tests whether a normal
# conditional distribution is adequate once the volatility clustering is
# modelled. Prices download at run time, so there is no data file to fetch.
#
#   Rscript R/02_sp500_arma_garch.R

library(quantmod)
library(tseries)
library(randtests)
library(fGarch)

set.seed(123)

dir.create("output", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

START <- "2000-01-01"
END   <- "2022-11-12"

log_lines <- character(0)
say <- function(...) {
  line <- paste0(...)
  cat(line, "\n", sep = "")
  log_lines <<- c(log_lines, line)
}

# ---- 1. data ---------------------------------------------------------------

getSymbols("^GSPC", from = START, to = END, auto.assign = TRUE)
returns <- na.omit(as.numeric(dailyReturn(GSPC, type = "log")))

say("S&P 500 daily log returns, ", START, " to ", END)
say("  observations: ", length(returns))
say("  mean:         ", sprintf("%.6f", mean(returns)))
say("  sd:           ", sprintf("%.6f", sd(returns)))

png("figures/log_returns.png", width = 1000, height = 500)
plot(returns, type = "l", col = "steelblue",
     main = "S&P 500 daily log returns", xlab = "Observation", ylab = "Log return")
abline(h = 0, col = "red", lty = 2)
dev.off()

# ---- 2. dependence and stationarity ----------------------------------------

runs <- runs.test(returns)
say("")
say("Runs test:            p = ", sprintf("%.4g", runs$p.value))

lb <- Box.test(returns, lag = 10, type = "Ljung-Box")
say("Ljung-Box (lag 10):   p = ", sprintf("%.4g", lb$p.value))

adf <- suppressWarnings(adf.test(returns, alternative = "stationary"))
say("ADF:                  stat = ", sprintf("%.3f", adf$statistic),
    ", p = ", sprintf("%.4g", adf$p.value))

# ---- 3. ARMA order by AIC --------------------------------------------------

arma_aic <- matrix(NA_real_, nrow = 4, ncol = 4,
                   dimnames = list(paste0("p", 0:3), paste0("q", 0:3)))
for (p in 0:3) {
  for (q in 0:3) {
    fit <- try(arima(returns, order = c(p, 0, q)), silent = TRUE)
    if (!inherits(fit, "try-error")) arma_aic[p + 1, q + 1] <- AIC(fit)
  }
}
best <- which(arma_aic == min(arma_aic, na.rm = TRUE), arr.ind = TRUE)[1, ]
arma_p <- best[["row"]] - 1
arma_q <- best[["col"]] - 1

say("")
say("Best ARMA by AIC:     ARMA(", arma_p, ",", arma_q, ")  AIC = ",
    sprintf("%.2f", min(arma_aic, na.rm = TRUE)))

# ---- 4. ARMA-GARCH order by AIC --------------------------------------------
# fGarch divides its information criteria by n, so these are comparable with
# each other but not with the arima() AIC above.

fits <- list()
aics <- c()
for (i in 0:2) for (j in 0:2) for (p in 1:2) for (q in 0:2) {
  tag <- sprintf("arma(%d,%d)+garch(%d,%d)", i, j, p, q)
  fit <- try(
    garchFit(substitute(~ arma(i, j) + garch(p, q),
                        list(i = i, j = j, p = p, q = q)),
             data = returns, trace = FALSE),
    silent = TRUE
  )
  if (!inherits(fit, "try-error")) {
    fits[[tag]] <- fit
    aics[tag]   <- fit@fit$ics[["AIC"]]
  }
}

best_tag   <- names(which.min(aics))
best_model <- fits[[best_tag]]

say("")
say("Models fitted:        ", length(fits), " of 54")
say("Best ARMA-GARCH:      ", best_tag, "  AIC = ", sprintf("%.5f", aics[[best_tag]]))
say("Runner up:            ", names(sort(aics)[2]), "  AIC = ",
    sprintf("%.5f", sort(aics)[2]))

# ---- 5. residual diagnostics -----------------------------------------------

resid_std <- na.omit(residuals(best_model, standardize = TRUE))

png("figures/qq_normal.png", width = 600, height = 600)
qqnorm(resid_std, main = "Standardised residuals against the normal")
qqline(resid_std, col = "red")
dev.off()

lb_resid <- Box.test(resid_std, lag = 10, type = "Ljung-Box")
say("")
say("Ljung-Box on standardised residuals (lag 10): p = ",
    sprintf("%.4g", lb_resid$p.value))
say("Excess kurtosis of standardised residuals:    ",
    sprintf("%.3f", mean((resid_std - mean(resid_std))^4) / sd(resid_std)^4 - 3))

# ---- 6. skewed generalised error distribution ------------------------------

sged_model <- try(
  garchFit(~ arma(1, 1) + garch(2, 2), data = returns,
           cond.dist = "sged", trace = FALSE),
  silent = TRUE
)

if (!inherits(sged_model, "try-error")) {
  say("")
  say("ARMA(1,1)+GARCH(2,2) under SGED")
  say("  AIC, normal: ", sprintf("%.5f", fits[["arma(1,1)+garch(2,2)"]]@fit$ics[["AIC"]]))
  say("  AIC, SGED:   ", sprintf("%.5f", sged_model@fit$ics[["AIC"]]))
  say("  shape:       ", sprintf("%.4f", sged_model@fit$par[["shape"]]))
  say("  skew:        ", sprintf("%.4f", sged_model@fit$par[["skew"]]))

  resid_sged <- na.omit(residuals(sged_model, standardize = TRUE))
  png("figures/qq_sged.png", width = 600, height = 600)
  qqplot(qsged(ppoints(length(resid_sged)),
               nu = sged_model@fit$par[["shape"]],
               xi = sged_model@fit$par[["skew"]]),
         sort(resid_sged),
         main = "Standardised residuals against the SGED",
         xlab = "Theoretical quantiles", ylab = "Sample quantiles")
  abline(0, 1, col = "red")
  dev.off()
} else {
  say("")
  say("SGED fit failed to converge.")
}

writeLines(log_lines, "output/results.txt")
