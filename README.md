# time-series-econometrics

ARMA-GARCH models of S&P 500 daily returns in R. It started as a homework for Statistical Methods in Finance (IE 522) at the University of Illinois. I later rewrote it and added a test for asymmetry, out-of-sample forecasts and a value-at-risk backtest.

## Data

Daily log returns from January 2000 to November 2022, 5,754 observations, downloaded from Yahoo Finance when a script runs. They are not independent. A runs test rejects at p = 3.5e-06 and Ljung-Box at ten lags rejects at p < 1e-15. The ADF statistic is -18.1, so the series can be modelled in levels.

## Model selection

`R/02_sp500_arma_garch.R`. For the mean alone, AIC picks ARMA(0,1). With a GARCH term added, I searched ARMA(0-2, 0-2) with GARCH(1-2, 0-2) and all 54 models converged. ARMA(1,1)-GARCH(2,2) has the lowest AIC at -6.43447, and ARMA(0,2)-GARCH(2,2) is at -6.43445. Two parts in 100,000 is not a difference, so I don't read anything into the choice between them.

The standardised residuals pass Ljung-Box (p = 0.53), so the serial dependence is gone. Their excess kurtosis is still 1.73. Refitting under a skewed generalised error distribution improves the AIC to -6.49688, with a shape of 1.37 against 2 for the normal and a skew of 0.86, which means a longer left tail.

## Asymmetry

`R/03_gjr_garch.R`. GJR-GARCH(1,1) adds a term that only switches on after a negative return. Under the same distribution, its coefficient is 0.169 (p = 0.004) and the likelihood ratio against plain GARCH(1,1) is 163.5 on one degree of freedom. A fall raises next-day variance much more than a rise of the same size, which `figures/news_impact.png` shows.

## Out-of-sample forecasts

`R/04_rolling_forecast.R`. Since AIC could not separate the top two models, I checked whether either forecasts better. Starting at observation 2,001, each model predicts the next day's variance on an expanding window, refitted every 250 days, which gives 3,754 forecasts scored against the squared return.

| Model | MSE | QLIKE |
|---|---|---|
| ARMA(1,1)-GARCH(2,2) | 5.312e-07 | -6.4625 |
| ARMA(0,2)-GARCH(2,2) | 5.303e-07 | -6.4422 |
| GARCH(1,1), no ARMA term | 5.466e-07 | -6.4646 |

Lower is better on both. The two AIC winners split, one ahead on MSE and the other on QLIKE, and the plain GARCH(1,1) has the best QLIKE of the three. QLIKE is the loss that is robust to a noisy variance proxy like the squared return, and on that loss the extra lags do not improve the forecast.

## Value at risk

`R/05_var_backtest.R`. One-day VaR from ARMA(1,1)-GARCH(2,2) under each distribution, scored with Kupiec's coverage test and Christoffersen's independence test.

Under the normal, losses beat the 1% VaR on 134 of 5,754 days, 2.3%, and Kupiec rejects outright. Under SGED it is 69 days, 1.2%, and Kupiec does not reject (p = 0.14). Christoffersen still rejects at p = 0.010, so the breaches arrive in clusters. At 5% the SGED model passes both tests, with 5.3% of days and p-values of 0.30 and 0.053.

Both models are fitted on the full sample and backtested on it, so these numbers flatter both of them. The comparison between the two distributions is still fair, since each gets the same advantage.

## Running it

```
Rscript R/01_garch_simulation.R
Rscript R/02_sp500_arma_garch.R
Rscript R/03_gjr_garch.R
Rscript R/04_rolling_forecast.R
Rscript R/05_var_backtest.R
```

They need `quantmod`, `tseries`, `randtests`, `fGarch` and `rugarch`. Each script writes its results to `output/` and its plots to `figures/`.

`01_garch_simulation.R` checks the simulator rather than the data. A GARCH(2,2) with persistence 0.95 implies an unconditional variance of 2.0, and 10,000 simulated draws give 2.035.

`fGarch` reports information criteria already divided by the sample size, so the GARCH AIC values above can be compared with each other but not with the ARMA AIC from `arima()`.
