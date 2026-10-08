# Asymmetry in the variance equation. The SGED skew came out below one, so this
# checks whether negative shocks also raise conditional variance more than
# positive ones. Both fits use rugarch.

library(quantmod)
library(rugarch)

set.seed(123)

dir.create("output", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
sink("output/gjr_garch.txt", split = TRUE)

invisible(getSymbols(Symbols = "^GSPC", from = "2000-01-01", to = "2022-11-12", auto.assign = TRUE))
log_returns <- na.omit(as.numeric(dailyReturn(GSPC, type = "log")))

garch_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(1, 1), include.mean = TRUE),
  distribution.model = "sged"
)

# gamma is the asymmetry term
gjr_spec <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(1, 1), include.mean = TRUE),
  distribution.model = "sged"
)

garch_fit <- ugarchfit(garch_spec, data = log_returns, solver = "hybrid")
gjr_fit <- ugarchfit(gjr_spec, data = log_returns, solver = "hybrid")

garch_ic <- infocriteria(garch_fit)
gjr_ic <- infocriteria(gjr_fit)

cat("GARCH(1,1) SGED\n")
cat("  AIC:", round(garch_ic["Akaike", ], 5), "\n")
cat("  BIC:", round(garch_ic["Bayes", ], 5), "\n")
cat("  loglik:", round(likelihood(garch_fit), 2), "\n\n")

cat("GJR-GARCH(1,1) SGED\n")
cat("  AIC:", round(gjr_ic["Akaike", ], 5), "\n")
cat("  BIC:", round(gjr_ic["Bayes", ], 5), "\n")
cat("  loglik:", round(likelihood(gjr_fit), 2), "\n\n")

coefs <- gjr_fit@fit$matcoef
print(round(coefs[grep("gamma", rownames(coefs)), , drop = FALSE], 5))

# GJR nests GARCH at gamma = 0, so this is an LR test with 1 df
lr <- 2 * (likelihood(gjr_fit) - likelihood(garch_fit))
cat("\nLR, GJR vs GARCH:", round(lr, 3), "p:", 1 - pchisq(lr, df = 1), "\n")

# News impact curve, symmetric under GARCH and kinked under GJR
ni1 <- newsimpact(garch_fit)
ni2 <- newsimpact(gjr_fit)

png("figures/news_impact.png", width = 800, height = 500)
plot(ni1$zx, ni1$zy, type = "l", col = "blue",
     ylim = range(c(ni1$zy, ni2$zy)),
     main = "News impact curve",
     xlab = "Shock", ylab = "Conditional variance")
lines(ni2$zx, ni2$zy, col = "red")
legend("top", legend = c("GARCH(1,1)", "GJR-GARCH(1,1)"),
       col = c("blue", "red"), lty = 1, bty = "n")
dev.off()

sink()
