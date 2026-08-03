# Threshold-VAR
Example implementation of a Bayesian threshold VAR for studying nonlinear, state-dependent macroeconomic effects using regime-specific and evolving-state impulse responses.

The model allows the VAR coefficients and shock covariance matrix to differ across low- and high-capacity-utilization regimes and considers both fixed-state and evolving-state generalised impulse responses.

## Repository Structure

- **`Data`** — Quarterly macroeconomic data and transformed series used in the estimation.

- **`Code`** — MATLAB scripts for data preparation, regime classification, prior specification, Gibbs sampling, and impulse-response computation.

- **`Model`** — Description of the threshold VAR specification, Bayesian estimation procedure, prior assumptions, and the construction of fixed-state and evolving-state impulse responses.
