# Simulates a GARCH(2,2) with known parameters and compares the sample moments
# with the ones the model implies. Sanity check on the simulator.

library(rugarch)

omega  <- 0.1
alpha1 <- 0.15
alpha2 <- 0.10
beta1  <- 0.50
beta2  <- 0.20
n      <- 10000

persistence <- alpha1 + alpha2 + beta1 + beta2
stopifnot(persistence < 1)   # stationarity

spec <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(2, 2)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm",
  fixed.pars = list(omega = omega,
                    alpha1 = alpha1, alpha2 = alpha2,
                    beta1  = beta1,  beta2  = beta2)
)

sim     <- ugarchpath(spec, n.sim = n, rseed = 123)
a_t     <- as.numeric(fitted(sim))
sigma_t <- as.numeric(sigma(sim))

# Mean is zero by construction: no constant in the mean equation and zero-mean
# innovations. Unconditional variance is omega / (1 - sum(alpha) - sum(beta)).
implied <- omega / (1 - persistence)

out <- c(
  sprintf("GARCH(2,2), %d observations", n),
  sprintf("persistence                    %.4f", persistence),
  sprintf("sample mean of a_t             %.5f", mean(a_t)),
  sprintf("implied unconditional variance %.4f", implied),
  sprintf("sample variance of a_t         %.4f", var(a_t)),
  sprintf("mean conditional variance      %.4f", mean(sigma_t^2))
)

writeLines(out)
dir.create("output", showWarnings = FALSE)
writeLines(out, "output/garch_simulation.txt")
