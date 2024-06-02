The GARCH folder contains the commands and functions needed to run a multivariate constant correlation GARCH(1,1) SSM model. The data generation process, generate a 4-d model, which corresponds to a similar analysis done in the manuscript, however, this can be changed as intended.

Brief description of the functions.

MvtKalFiltGarch
Details: Computes the Kalman-filter from a state space model with a constant correlation GARCH(1,1) model. <br>
Arguments: Data, initializing constants for the state means and variances, and observational level variance and GARCH parameters. <br>
Returns: Filtered state means and variance, 1-step ahead forecasts and smoothed means and variances.

MvtKalFFBS
Details: Computes the backward-sampled states means and variances obtained once a Kalman-filter Iteration has run.
Arguments: Output from a Kalman Filter.
Returns: Backward-sampled states means and variances.

MvtLkhd
Details: Computes the log likelihood of the data given the parameter values.
Arguments: Data, GARCH parameters.
Returns: Log of the likelihood of the model


dmvnormC
Details: Fast Rcpp computation of the multivariate normal.
Arguments: Data, means, covariances.
Output: log of the density function.

MVT_Garch_waic1
Details: Computes the WAIC based on the posterior estimates.
Arguments: Data, posterior estimates.
Output: WAIC.
