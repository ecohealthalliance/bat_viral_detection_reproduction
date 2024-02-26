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
  
  vector[3] beta; // three betas, representing community-level intercept, pregnancy effect, and lactation effect
  
   // parameters related to varying effects
    
  cholesky_factor_corr[3] L_Omega; // Cholesky correlation matrix
  vector<lower=0>[3] sigma_host_species; // sd for the species-specific int/preg/lac effects
  matrix[3, N_host_species] z; // matrix of standardized int/preg/lac offsets
  
  vector[N_year] alpha_year_tilde;
  vector[N_country] alpha_country_tilde;
  vector[N_specimen_type_group] alpha_specimen_type_group_tilde;
  vector[N_diagnostic_laboratory_name] alpha_diagnostic_laboratory_name_tilde;
  
  vector<lower=0>[4] sigma_vector;
}

///////////////////////////////////////////////////////////////////////////////
  
  
transformed parameters {
  
  // create matrix of int/preg/lac offsets
  
  matrix[3, N_host_species] z_s; 
  z_s = diag_pre_multiply(sigma_host_species, L_Omega) * z;
  
  // create matrix of species-specific int/preg/lac effects
  
  matrix[3, N_host_species] beta_host_species;
  for (i in 1:N_host_species) { // for each host species:
    
    beta_host_species[1, i] = beta[1] + z_s[1, i]; // build intercept
    beta_host_species[2, i] = beta[2] + z_s[2, i]; // build pregnancy effect
    beta_host_species[3, i] = beta[3] + z_s[3, i]; // build lactation effect
  }
  
  // transform non-centered varying effects
  
  vector[N_year] alpha_year;
  vector[N_country] alpha_country;
  vector[N_specimen_type_group] alpha_specimen_type_group;
  vector[N_diagnostic_laboratory_name] alpha_diagnostic_laboratory_name;
  
  for (i in 1:N_year) {
   
    alpha_year[i] = 
    sigma_vector[1] * alpha_year_tilde[i];
  }
  
  for (i in 1:N_country) {
    
    alpha_country[i] = 
    sigma_vector[2] * alpha_country_tilde[i];
  }
    
  for (i in 1:N_specimen_type_group) {
   
    alpha_specimen_type_group[i] = 
    sigma_vector[3] * alpha_specimen_type_group_tilde[i];
  }
    
  for (i in 1:N_diagnostic_laboratory_name) {
    
    alpha_diagnostic_laboratory_name[i] = 
    sigma_vector[4] * alpha_diagnostic_laboratory_name_tilde[i];
  }
  
  // calculate overall linear predictor
  
  vector[N] alpha;
  
  for (i in 1:N) {
      
    alpha[i] = 
    
    // species-specific intercept
    beta_host_species[1, host_species[i]] + 
    
    // species-specific pregnancy effect
    beta_host_species[2, host_species[i]] * pregnant_mod[i] + 
    
    // species-specific lactation effect
    beta_host_species[3, host_species[i]] * lactating_mod[i] + 
    
    // varying intercept values for all other clusters
    alpha_year[year[i]] + 
    alpha_country[country[i]] + 
    alpha_specimen_type_group[specimen_type_group[i]] +
    alpha_diagnostic_laboratory_name[diagnostic_laboratory_name[i]];
  }
}

///////////////////////////////////////////////////////////////////////////////
  
  
model {
    

  // priors for main effects
  
  beta[1] ~ normal(0, 2);
  beta[2] ~ normal(0, 1);
  beta[3] ~ normal(0, 1);
  
  // priors for varying effects
    
  L_Omega ~ lkj_corr_cholesky(2);
  sigma_host_species ~ exponential(1);
  to_vector(z) ~ normal(0, 1);
  
  alpha_year_tilde ~ normal(0, 1);
  alpha_country_tilde ~ normal(0, 1);
  alpha_specimen_type_group_tilde ~ normal(0, 1);
  alpha_diagnostic_laboratory_name_tilde ~ normal(0, 1);
   
  sigma_vector ~ exponential(1);
  
  // likelihood
    
  virus_detected ~ bernoulli_logit(alpha);
}

///////////////////////////////////////////////////////////////////////////////


generated quantities {
  
  // output correlation matrix for intercept and slope parameters
  
  matrix[3, 3] Omega;
  Omega = multiply_lower_tri_self_transpose(L_Omega);
}
