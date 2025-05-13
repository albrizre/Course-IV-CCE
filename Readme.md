**Materials** for the course *Bayesian modeling of spatio-temporal point pattern data* taught during the IV Workshop of the Colombian Statistical Society in May 2025.

The *main folder* contains the codes to run the different models described during the course, including data loading and preparation, definition of constants and initial values, and model fitting with **nimble**. In particular:

- "0 Data exploration.R" contains some lines of code for an initial exploration of the data.
- "1 Homogeneous Poisson model.R" contains the code for fitting the Poisson homogeneous model.
- "2 Inhomogeneous Poisson model v1.R", "2 Inhomogeneous Poisson model v2.R", and "2 Inhomogeneous Poisson model v3.R" contain the code for fitting the following three versions of an inhomogeneous Poisson model: including spatial covariates (v1), including spatial covariates and a temporal effect decaying from day t = 200 (v2), and including spatial covariates and a first-order random walk for the temporal effect (v3).
- "3 Splines for the spatial intensity.R", "3 Underdetection events.R", "3 Temporal location uncertainty.R" contain the code for fitting the "advanced" examples described in the course.
-  "4 Residual analysis.R" contains the code for doing a spatial residual analysis on the Poisson homogeneous model (note that the output from this model needs to be in the "Outputs" folder, so run code "1 Homogeneous Poisson model.R" or any other of the codes to perform the residual analysis).

**Furthermore:**

- Folder *Models* contains the **nimble** codes of each of the Bayesian point process models explained.
- Folder *Data* contains the data required to use the codes. Note that to avoid violating the privacy of the data, the time locations have been obtained from a uniform distribution on [0,365] (the results shown in the slides correspond to the real dataset).
- Folder *Outputs* should be used for saving the outputs of the fitted models according to the codes. The folder is presented as empty because the outputs are quite large. Be careful with the number of interations/chains you choose for fitting the model.
- Folder *Slides* contains both the slides of the course and the slides of the talk given during the workshop about an extension of the self-exciting model.
