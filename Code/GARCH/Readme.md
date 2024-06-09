The GARCH folder contains the commands and functions needed to run a multivariate constant correlation GARCH(1,1) SSM model. The data generation process generates a 4-d model, which corresponds to a similar analysis done in the manuscript, however, this can be changed as intended. 

The R-script, step_by_step_guide.R, guides the user from model generation to estimation and residual analysis.

Brief description of some of the key functions.

**MvtKalFiltGarch**
Details: Computes the Kalman-filter from a state space model with a constant correlation GARCH(1,1) model. <br>
Arguments: Data, initializing constants for the state means and variances, and observational level variance and GARCH parameters. <br>
Returns: Filtered state means and variance, 1-step ahead forecasts and smoothed means and variances.  <br>

**MvtKalFFBS**
Details: Computes the backward-sampled states means and variances obtained once a Kalman-filter Iteration has run.  <br>
Arguments: Output from a Kalman Filter.  <br>
Returns: Backward-sampled states means and variances.  <br>

**MvtLkhd**  <br>
Details: Computes the log likelihood of the data given the parameter values.  <br>
Arguments: Data, GARCH parameters.  <br>
Returns: Log of the likelihood of the model.  <br>

**MVT_Garch_waic1**  <br>
Details: Computes the WAIC based on the posterior estimates.  <br>
Arguments: Data, posterior estimates.  <br>
Output: WAIC.  <br>
