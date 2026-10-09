# Simulates a GARCH(2,2) with known parameters and checks that the sample
# variance is close to the variance the model implies. Fails with a non-zero
# exit if it is not.
#
#   Rscript tests/check_simulation.R

suppressPackageStartupMessages(library(rugarch))

pars <- list(omega = 0.1, alpha1 = 0.15, alpha2 = 0.10, beta1 = 0.50, beta2 = 0.20)
spec <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(2, 2)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm",
  fixed.pars         = pars
)
sim <- ugarchpath(spec, n.sim = 20000, rseed = 7)
a_t <- as.numeric(fitted(sim))

implied <- pars$omega / (1 - pars$alpha1 - pars$alpha2 - pars$beta1 - pars$beta2)
ratio <- var(a_t) / implied
cat(sprintf("implied %.3f  sample %.3f  ratio %.3f\n", implied, var(a_t), ratio))
if (abs(ratio - 1) > 0.1 || abs(mean(a_t)) > 0.1) {
  stop("simulated moments do not match the model")
}
cat("pass\n")
