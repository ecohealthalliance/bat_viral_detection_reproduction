This repository contains code, data, and figure/table files for an analysis of reproductive effects on bat viral detection that leverages data from the [USAID-funded PREDICT project](https://ohi.vetmed.ucdavis.edu/programs-projects/predict-project). More specifically, data for this project were accessed from the [Emering Infectious Disease Information Technology Hub](https://www.eidith.org/) (EIDITH) database using the [`eidith` R package](https://ecohealthalliance.github.io/eidith/).

### Repository contents

The contents of the repository subdirectories are as follows:

- `R/` contains internal functions used in the repository code

- `data/` contains various lookup tables used during data cleaning as well as the cleaned data for model-fitting

- `outputs/` contains all figures and tables associated with the analysis

- `scripts/` contains the core repository analyses:
  - `01_data_prep.R` accesses, cleans, and filters EIDITH data for further analysis
  - `02_stan_model_fitting.R` fits hierarchical Bayesian models to the cleaned datasets using `rstan`
  - `03_plotting.R` generates all figures

- `stan/` contains cleaned data subsets in Stan-friendly formats and the Stan models used in the analysis

### Reproducing the analysis

In order to facilitate reproducibility for those without EIDITH access, cleaned data analyzed here are contained in both the `data/cleaned_data` subdirectory (full analysis set in CSV format) and the `stan/cleaned_data` subdirectory (full analysis set and data subsets as Stan-friendly list objects). Users cloning the repository will therefore be able to proceed with the recreation of the analysis starting with `scripts/02_stan_model_fitting.R`, which will fit the Stan models, saving them to `stan/saved_models`. Subsequently, `03_plotting.R` will recreate all figures using this saved model output.