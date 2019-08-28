data {
  
  
  // number of observations
  
  int<lower=1> N;
  
  
  // outcome variable
  
  int<lower=0, upper=1> virus_detected[N];
  
  
  // variables related to varying effects
  
  int<lower=1> N_host_species;
  int<lower=1, upper=N_host_species> host_species[N];
  
  int<lower=1> N_year;
  int<lower=1, upper=N_year> year[N];
  
  int<lower=1> N_country;
  int<lower=1, upper=N_country> country[N];
  
  int<lower=1> N_specimen_type_group;
  int<lower=1, upper=N_specimen_type_group> specimen_type_group[N];
  
  int<lower=1> N_diagnostic_laboratory_name;
  int<lower=1, upper=N_diagnostic_laboratory_name> diagnostic_laboratory_name[N];
  
  
  // variables related to reproductive effects

  int<lower=0, upper=1> pregnant_mod[N];
  
  int<lower=0, upper=1> lactating_mod[N];
}

///////////////////////////////////////////////////////////////////////////////
  
  
parameters {
    
    
  // parameters related to main effects
  
  real mu_alpha;
  
  real beta_pregnant_mod;
  
  real beta_lactating_mod;
  
  
  // parameters related to varying effects
    
  vector[N_host_species] alpha_host_species_offset_tilde;
  vector[N_year] alpha_year_tilde;
  vector[N_country] alpha_country_tilde;
  vector[N_specimen_type_group] alpha_specimen_type_group_tilde;
  vector[N_diagnostic_laboratory_name] alpha_diagnostic_laboratory_name_tilde;
  
  vector<lower=0>[5] sigma_vector;
  real<lower=0> scale_for_sigmas;
}

///////////////////////////////////////////////////////////////////////////////
  
  
transformed parameters {
  
  
  vector[N_host_species] alpha_host_species_offset;
  vector[N_year] alpha_year;
  vector[N_country] alpha_country;
  vector[N_specimen_type_group] alpha_specimen_type_group;
  vector[N_diagnostic_laboratory_name] alpha_diagnostic_laboratory_name;
  
  vector[N] alpha;
  
  
  // transform non-centered varying effects
  
  for (i in 1:N_host_species) {
   
    alpha_host_species_offset[i] = 
    sigma_vector[1] * alpha_host_species_offset_tilde[i];
  }
  
  for (i in 1:N_year) {
   
    alpha_year[i] = 
    sigma_vector[2] * alpha_year_tilde[i];
  }
  
  for (i in 1:N_country) {
    
    alpha_country[i] = 
    sigma_vector[3] * alpha_country_tilde[i];
  }
    
  for (i in 1:N_specimen_type_group) {
   
    alpha_specimen_type_group[i] = 
    sigma_vector[4] * alpha_specimen_type_group_tilde[i];
  }
    
  for (i in 1:N_diagnostic_laboratory_name) {
    
    alpha_diagnostic_laboratory_name[i] = 
    sigma_vector[5] * alpha_diagnostic_laboratory_name_tilde[i];
  }
    
  
  // calculate overall linear predictor
  
  for (i in 1:N) {
      
    alpha[i] = 
    
    mu_alpha + 
    
    alpha_host_species_offset[host_species[i]] +
    alpha_year[year[i]] + 
    alpha_country[country[i]] + 
    alpha_specimen_type_group[specimen_type_group[i]] +
    alpha_diagnostic_laboratory_name[diagnostic_laboratory_name[i]] +
    
    (beta_pregnant_mod * pregnant_mod[i]) +
    (beta_lactating_mod * lactating_mod[i]);
  }
}

///////////////////////////////////////////////////////////////////////////////
  
  
model {
    

  // priors for main effects
  
  mu_alpha ~ normal(0, 5);
  
  beta_pregnant_mod ~ normal(0, 5);
  
  beta_lactating_mod ~ normal(0, 5);
  
  
  // priors for varying effects
    
  alpha_host_species_offset_tilde ~ normal(0, 1);
  alpha_year_tilde ~ normal(0, 1);
  alpha_country_tilde ~ normal(0, 1);
  alpha_specimen_type_group_tilde ~ normal(0, 1);
  alpha_diagnostic_laboratory_name_tilde ~ normal(0, 1);
   
  sigma_vector ~ cauchy(0, scale_for_sigmas);
  scale_for_sigmas ~ exponential(1);
  
    
  // likelihood
    
  virus_detected ~ bernoulli_logit(alpha);
}

///////////////////////////////////////////////////////////////////////////////


generated quantities {
  
  
  vector[N] log_lik;
  
  vector[N_host_species] alpha_host_species;
  
  
  for (i in 1:N) {

    log_lik[i] = bernoulli_logit_lpmf(virus_detected[i] | alpha[i]);
  }
  
  
  for (i in 1:N_host_species) {
   
    alpha_host_species[i] = mu_alpha + alpha_host_species_offset[i];
  }
}
