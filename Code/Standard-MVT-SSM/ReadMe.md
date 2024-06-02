The Standard-MVT-SSM folder contains the commands and functions needed to run a standard SSM model.
In the data generation process, we use generate independent series, however, this can be modified as needed.


Bried description of the functions.


**MvtKalFilt**

Details: Computes the Kalman-filter from a state space model. <be>
Arguments: Data, initializing constants for the state means and variances, and observational level variance.
Returns: Filtered state means and variance, 1-step ahead forecasts and smoothed means and variances.


**MvtKalFFBS**

Details: Computes the backward-sampled states means and variances obtained once a Kalman-filter Iteration has run.
Arguments: Output from a Kalman Filter.
Returns: Backward-sampled states means and variances.

**dmvnormC**

Details: Fast Rcpp computation of the multivariate normal.
Arguments: Data, means, covariances.
Output: log of the density function.

**waic1_std**

Details: Computes the WAIC based on the posterior estimates.
Arguments: Data, posterior estimates.
Output: WAIC.
