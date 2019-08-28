#' ---
#' title: "Bat Viral Detection and Reproduction Analysis - Stan Model Fitting"
#' author: "EcoHealth Alliance M&A Team - Code Drafted by Evan Eskew"
#' ---

# /*
#==============================================================================
# */


#+ loading_chunk, echo=FALSE, message=FALSE, results="hide"


# Load packages

library(rprojroot)
library(rstan)
library(dplyr)
library(stringr)


# Set working directory

setwd(find_rstudio_root_file())


# Get data files

data.files <- list.files("stan/cleaned_data/")


# Set Stan model parameters

iter <- 3500
warmup <- 1000
chains <- 4
cores <- 4
adapt_delta <- 0.99
stepsize <- 0.5


# Fit Stan models

dataset.model.link <-
  c("stan/model_code/bat_reproduction.stan",
    "stan/model_code/bat_reproduction.stan",
    "stan/model_code/bat_reproduction.stan",
    "stan/model_code/bat_reproduction.stan",
    "stan/model_code/bat_reproduction_no_test_variation.stan",
    "stan/model_code/bat_reproduction.stan"
  )

assert_that(length(data.files) == length(dataset.model.link))

for(i in 1:length(data.files)) {
  
  # Assign one dataset
  
  stan.data <- readRDS(paste0("stan/cleaned_data/", data.files[i]))
  
  # Fit the model 
  
  set.seed(8)
  
  fit.model <- stan(
    file = dataset.model.link[i], 
    data = stan.data, iter = iter, warmup = warmup,
    chains = chains, cores = cores, verbose = TRUE,
    control = list(adapt_delta = adapt_delta, stepsize = stepsize)
  )
  
  # Save the model
  
  saveRDS(
    fit.model,
    file = paste0(
      "stan/saved_models/",
      data.files[i] %>%
        str_replace("dat.", "model.") %>%
        str_replace(".stan.rds", ""),
      ".rds")
  )
}
