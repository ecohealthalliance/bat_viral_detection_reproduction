get_stan_data <- function(x) {
  
  return(
    
    list(
      
      N = nrow(x),
      
      virus_detected = x$virus_detected,
      
      lactating_mod = to_stan_factor(x$lactating_mod) - 1,
      
      pregnant_mod = to_stan_factor(x$pregnant_mod) - 1,
      
      year = to_stan_factor(x$year),
      N_year = stan_factor_count(x$year),
      
      country = to_stan_factor(x$country),
      N_country = stan_factor_count(x$country),
      
      specimen_type_group = to_stan_factor(x$specimen_type_group),
      N_specimen_type_group = stan_factor_count(x$specimen_type_group),
      
      diagnostic_laboratory_name = 
        to_stan_factor(x$diagnostic_laboratory_name),
      N_diagnostic_laboratory_name = 
        stan_factor_count(x$diagnostic_laboratory_name),
      
      test_requested = to_stan_factor(x$test_requested),
      N_test_requested = stan_factor_count(x$test_requested),
      
      test_requested_mod = to_stan_factor(x$test_requested_mod),
      N_test_requested_mod = stan_factor_count(x$test_requested_mod),
      
      host_order = to_stan_factor(x$order),
      N_host_order = stan_factor_count(x$order),
      
      host_family = to_stan_factor(x$family),
      N_host_family = stan_factor_count(x$family),
      
      host_genus = to_stan_factor(x$genus),
      N_host_genus = stan_factor_count(x$genus),
      
      host_species = to_stan_factor(x$binomial),
      N_host_species = stan_factor_count(x$binomial)
    )
  )
}
