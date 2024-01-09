#' ---
#' title: "Bat Viral Detection and Reproduction Analysis - Stan Model Fitting"
#' author: "EcoHealth Alliance M&A Team - Code Drafted by Evan Eskew"
#' ---

# /*
#==============================================================================
# */


#+ loading_chunk, echo=FALSE, message=FALSE, results="hide"


# Loop through all datasets twice, the first time fitting the pooled effect
# model, the second time fitting the varying intercepts and slopes model


# Load packages

library(cmdstanr)
library(dplyr)
library(stringr)
library(assertthat)


# Get data files

data.files <- rep(list.files("stan/cleaned_data/"), times = 2)


# Set Stan model parameters

chains <- 4
iter_warmup <- 1000
iter_sampling <- 2500
adapt_delta <- 0.99
stepsize <- 0.5


# Fit Stan models

dataset.model.link <-
  c(
    "stan/model_code/bat_reproduction.stan",
    "stan/model_code/bat_reproduction.stan",
    "stan/model_code/bat_reproduction.stan",
    "stan/model_code/bat_reproduction.stan",
    "stan/model_code/bat_reproduction_no_test_variation.stan",
    "stan/model_code/bat_reproduction.stan",
    "stan/model_code/varying_ints_slopes_centered.stan",
    "stan/model_code/varying_ints_slopes_centered.stan",
    "stan/model_code/varying_ints_slopes_centered.stan",
    "stan/model_code/varying_ints_slopes_centered.stan",
    "stan/model_code/varying_ints_slopes_centered_no_test_variation.stan",
    "stan/model_code/varying_ints_slopes_centered.stan"
  )

assert_that(length(data.files) == length(dataset.model.link))

seeds <- rep(c(1, 8, 8, 8, 1, 8), times = 2)

# Loop through datasets and fit the models

for(i in 1:length(data.files)) {
  
  # Assign one dataset
  
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
      ifelse(
        str_detect(dataset.model.link[i], "varying"),
        # names for varying effect models
        data.files[i] %>%
          str_replace("dat.f", "model.f.v") %>%
          str_replace(".stan", ""),
        # names for pooled effect models
        data.files[i] %>%
          str_replace("dat.f", "model.f") %>%
          str_replace(".stan", "")
      )
    )
  )
}
