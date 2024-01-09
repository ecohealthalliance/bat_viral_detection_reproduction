data {
  
  
  // number of observations
  
  int<lower=1> N;
  
  
  // outcome variable
  
  array[N] int<lower=0, upper=1> virus_detected;
  
  
  // variables related to varying effects
  
  int<lower=1> N_host_species;
  array[N] int<lower=1, upper=N_host_species> host_species;
  
  int<lower=1> N_year;
  array[N] int<lower=1, upper=N_year> year;
  
  int<lower=1> N_country;
  array[N] int<lower=1, upper=N_country> country;
  
  int<lower=1> N_specimen_type_group;
  array[N] int<lower=1, upper=N_specimen_type_group> specimen_type_group;
  
  int<lower=1> N_test_requested_mod;
  array[N] int<lower=1, upper=N_test_requested_mod> test_requested_mod;
  
  int<lower=1> N_diagnostic_laboratory_name;
  array[N] int<lower=1, upper=N_diagnostic_laboratory_name> diagnostic_laboratory_name;
  
  
  // variables related to reproductive effects

  array[N] int<lower=0, upper=1> pregnant_mod;
  
  array[N] int<lower=0, upper=1> lactating_mod;
}

///////////////////////////////////////////////////////////////////////////////
  
  
parameters {
    
    
  // parameters related to main effects
  
  vector[3] beta;
  
   // parameters related to varying effects
    
  corr_matrix[3] omega;
  vector<lower=0>[3] sigma_host_species;
  array[N_host_species] vector[3] beta_host_species;
}

///////////////////////////////////////////////////////////////////////////////
  
  
transformed parameters {
  
  
  vector[N] alpha;
  
  // calculate overall linear predictor
  
  for (i in 1:N) {
      
    alpha[i] = 
    
     // species-specific intercept
    beta_host_species[host_species[i], 1] +
    
    // species-specific pregnancy effect
    beta_host_species[host_species[i], 2] * pregnant_mod[i] + 
    
    // species-specific lactation effect
    beta_host_species[host_species[i], 3] * lactating_mod[i];
  }
}

///////////////////////////////////////////////////////////////////////////////
  
  
model {
    

  // priors for main effects
  
  beta[1] ~ normal(0, 2);
  beta[2] ~ normal(0, 1);
  beta[3] ~ normal(0, 1);
  
  // priors for varying effects
    
  omega ~ lkj_corr(2);
  sigma_host_species ~ exponential(1);
  beta_host_species ~ multi_normal(beta, quad_form_diag(omega, sigma_host_species));
  
  // likelihood
    
  virus_detected ~ bernoulli_logit(alpha);
}
