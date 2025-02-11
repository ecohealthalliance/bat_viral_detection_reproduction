#' ---
#' title: "Bat Viral Detection and Reproduction Analysis - Stan Model Fitting"
#' author: "EcoHealth Alliance M&A Team - Code Drafted by Evan Eskew"
#' ---

# /*
#==============================================================================
# */


#+ loading_chunk, echo=FALSE, message=FALSE, results="hide"


# Loop through all modeling datasets three times to:
# 1) fit main effects models
# 2) fit varying intercepts models
# 3) fit varying intercepts and slopes models


# Load packages

library(cmdstanr)
library(dplyr)
library(stringr)
library(assertthat)


# Get data files

data.files <- rep(list.files("stan/cleaned_data/"), times = 3)


# Set Stan model parameters

chains <- 4
iter_warmup <- 1000
iter_sampling <- 2500
adapt_delta <- 0.99
stepsize <- 0.5


# Fit Stan models

dataset.model.link <-
  c(
    "stan/model_code/main_effects.stan",
    "stan/model_code/main_effects.stan",
    "stan/model_code/main_effects.stan",
    "stan/model_code/varying_ints.stan",
    "stan/model_code/varying_ints.stan",
    "stan/model_code/varying_ints.stan",
    "stan/model_code/varying_ints_slopes_non_centered.stan",
    "stan/model_code/varying_ints_slopes_non_centered.stan",
    "stan/model_code/varying_ints_slopes_non_centered.stan"
  )

assert_that(length(data.files) == length(dataset.model.link))

seeds <- rep(c(8, 8, 8), times = 3)

# Loop through datasets and fit the models

for(i in 1:length(data.files)) {
  
  # Assign dataset
  
  stan.data <- readRDS(paste0("stan/cleaned_data/", data.files[i]))
  
  # Fit the model 
  
  model <- cmdstan_model(dataset.model.link[i])
  
  fit.model <- model$sample(
    data = stan.data, 
    chains = chains,
    parallel_chains = chains,
    iter_warmup = iter_warmup, 
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta, 
    step_size = stepsize,
    seed = seeds[i]
  )
  
  # Save the model
  
  fit.model$save_object(
    file = paste0(
      "stan/saved_models/",
      case_when(
        str_detect(dataset.model.link[i], "main_effects") == TRUE ~
          # names for main effects models
          data.files[i] %>%
          str_replace("dat.f", "model.f.m") %>%
          str_replace(".stan", ""),
        str_detect(dataset.model.link[i], "varying_ints_slopes") == TRUE ~
          # names for varying intercepts and slopes models
          data.files[i] %>%
          str_replace("dat.f", "model.f.vis") %>%
          str_replace(".stan", ""),
        # names for varying intercepts models
        str_detect(dataset.model.link[i], "varying_ints") == TRUE ~ 
          data.files[i] %>%
          str_replace("dat.f", "model.f.vi") %>%
          str_replace(".stan", "")
      )
    )
  )
}
