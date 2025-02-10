data {
  
  
  // number of observations
  
  int<lower=1> N;
  
  
  // outcome variable
  
  array[N] int<lower=0, upper=1> virus_detected;
  
  
  // variables related to reproductive effects

  array[N] int<lower=0, upper=1> pregnant_mod;
  
  array[N] int<lower=0, upper=1> lactating_mod;
}

///////////////////////////////////////////////////////////////////////////////
  
  
parameters {
    
    
  // parameters related to main effects
  
  real mu_alpha;
  
  real beta_pregnant_mod;
  
  real beta_lactating_mod;
}

///////////////////////////////////////////////////////////////////////////////
  
  
transformed parameters {
  
  
  vector[N] alpha;
    
  
  // calculate overall linear predictor
  
  for (i in 1:N) {
      
    alpha[i] = 
    
    mu_alpha + 
    
    (beta_pregnant_mod * pregnant_mod[i]) +
    (beta_lactating_mod * lactating_mod[i]);
  }
}

///////////////////////////////////////////////////////////////////////////////
  
  
model {
    

  // priors for main effects
  
  mu_alpha ~ normal(0, 2);
  
  beta_pregnant_mod ~ normal(0, 1);
  
  beta_lactating_mod ~ normal(0, 1);
  
    
  // likelihood
    
  virus_detected ~ bernoulli_logit(alpha);
}
