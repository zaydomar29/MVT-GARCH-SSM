The Standard-MVT-SSM folder contains the commands and functions needed to run a standard SSM model.
In the data generation process, we use generate independent series, however, this can be modified as needed.


Bried description of the functions.


**MvtKalFilt** <br>
Details: Computes the Kalman-filter from a state space model. <br>
Arguments: Data, initializing constants for the state means and variances, and observational level variance. <br>
Returns: Filtered state means and variance, 1-step ahead forecasts and smoothed means and variances.


**MvtKalFFBS** <br>
Details: Computes the backward-sampled states means and variances obtained once a Kalman-filter Iteration has run. <br>
Arguments: Output from a Kalman Filter. <br>
Returns: Backward-sampled states means and variances. <br>

**dmvnormC** <br>
Details: Fast Rcpp computation of the multivariate normal. <br>
Arguments: Data, means, covariances. <br>
Output: log of the density function. <br>

**waic1_std** <br>
Details: Computes the WAIC based on the posterior estimates. <br>
Arguments: Data, posterior estimates. <br>
Output: WAIC. <br>
