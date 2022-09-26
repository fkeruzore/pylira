# pylira

Python wrapper for the `LIRA` library.
`LIRA` (LInear Regression in Astronomy, Sereno 2016) is an `R` library that performs linear regression using Bayesian inference and Gibbs sampling.
From the official documentation:

> [LIRA] Performs Bayesian linear regression and forecasting in astronomy. The method accounts for heteroscedastic errors in both the independent and the dependent variables, intrinsic scatters (in both variables) and scatter correlation, time evolution of slopes, normalization, scatters, Malmquist and Eddington bias, upper limits and break of linearity. The posterior distribution of the regression parameters is sampled with a Gibbs method exploiting the JAGS library.

See official documentation (https://cran.r-project.org/web/packages/lira/lira.pdf) and accompanying paper (Sereno 2016, https://ui.adsabs.harvard.edu/abs/2016MNRAS.455.2149S/abstract)