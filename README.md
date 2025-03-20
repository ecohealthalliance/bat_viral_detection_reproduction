# Reproductive effects on viral detection in bats

This repository contains code, data, figures, and tables that support the manuscript: 

Eskew, E.A., K.J. Olival, J.A.K. Mazet, and P. Daszak. 2025. [A global-scale dataset of bat viral detection suggests that pregnancy reduces viral shedding](https://doi.org/10.1098/rspb.2024.2381). Proceedings of the Royal Society B 292: 20242381. 

The paper's analyses broadly address reproductive effects on bat viral detection and leverage data from the [USAID-funded PREDICT project](https://ohi.vetmed.ucdavis.edu/programs-projects/predict-project). Data were accessed from the [Emering Infectious Disease Information Technology Hub](https://www.eidith.org/) (EIDITH) database using the [`eidith` R package](https://ecohealthalliance.github.io/eidith/).

---

### Repository Contents

The contents of the repository subdirectories are as follows:

- [`R`](/R/) contains internal functions used in the repository code

- [`data`](/data/) contains various lookup tables used during data cleaning as well as the cleaned viral detection data used for model-fitting

- [`outputs`](/outputs/) contains all figures and tables associated with the analyses

- [`scripts`](/scripts/) contains the core analysis scripts:
  - [`01_data_prep.R`](/scripts/01_data_prep.R) accesses, cleans, and filters EIDITH data for downstream analysis
  - [`02_stan_model_fitting.R`](/scripts/02_stan_model_fitting.R) fits Bayesian models to the cleaned datasets using [`cmdstanr`](https://mc-stan.org/cmdstanr/)
  - [`03_plotting.R`](/scripts/03_plotting.R) generates all figures

- [`stan`](/stan/) contains cleaned data in Stan-friendly formats and the Stan models used in the analysis

---

### Reproducing the Analyses

To facilitate reproducibility for those without EIDITH access, cleaned data analyzed here are contained in both the [`data/cleaned_data`](/data/cleaned_data/) subdirectory (full analysis set in CSV format) and the [`stan/cleaned_data`](/stan/cleaned_data/) subdirectory (full analysis set and data subsets as Stan-friendly objects). Users cloning the repository will therefore be able to proceed with the recreation of the analyses starting with [`02_stan_model_fitting.R`](/scripts/02_stan_model_fitting.R), which will fit the Stan models, saving them to `stan/saved_models`. Subsequently, [`03_plotting.R`](/scripts/03_plotting.R) will recreate all figures using the saved models.