# Clean workspace ----
rm(list = ls())

# Load functions  ----
devtools::load_all()
# source("R/model_functions.R") # Alternatively run this to load functions if not downloaded as an RStudio project

# Load packages ----
library(tidyverse)
library(deSolve)

# Define parameters ----
n_hhsize <- 3 # 2 # 1
n_pop_size <- 1000000
n_inf   <- 1
r_beta  <- 0.25 
n_contacts_comm <- 10 
r_sigma <- 0.35
r_gamma <- 0.2
n_contacts_hh <- 4
r_tau   <- 0.35
r_omega <- 0 #1/100
r_vax_omega <- 0 #1/100
redux_vax <- 0.1
eff_vax <- (1.0 - redux_vax)

n_exp_states <- 3
n_inf_states <- 2

r_dx          <- 0.10 # detection rate
p_alpha_dx    <- 0.20 # decrease in infectiousness
r_death       <- 0.00 #0.001 
p_death_inf   <- 0.02
r_growth_rate <- 1.00 # 1.01 # population growth rate
r_birth       <- r_death*r_growth_rate

max_time <- 200

times <- seq(0, max_time, by = 1)

# NPI policies ----
v_npi_times  <- c(0, 15, 16, 70, max_time + 1)
# v_npi_levels <- c(1, 1, 0.2, 0.2, 0.2) # Uncomment to try other NPI 
v_npi_levels <- c(1, 1, 0.01, 0.01, 0.01)
# v_npi_levels <- c(1, 1, 0.6, 0.6, 0.6)
# v_npi_levels <- c(1, 1, 1, 1, 1)

fun_npi <- approxfun(x = v_npi_times, y = v_npi_levels, method = "linear")

# Vaccination policies ----
v_time_stop_vax <- c(0, 40, 41, max_time + 1)
v_duration  <- diff(c(0, v_time_stop_vax))
# v_cum_prop_time <- c(0, 0.80, 0, 0) # Uncomment to try a vaccination strategy
v_cum_prop_time <- c(0, 0, 0, 0)

## Change of vaccination policies over time ---
v_vax_rates <- cbind(daily_rate(cum_prop = v_cum_prop_time[2], # We need constant rate to get up to cumulative coverage
                                duration = v_duration[2]),       # We need constant rate to get up to cumulative coverage
                     daily_rate(cum_prop = v_cum_prop_time[2],
                                duration = v_duration[2]),
                     daily_rate(cum_prop = v_cum_prop_time[3],
                                duration = v_duration[3]),
                     daily_rate(cum_prop = v_cum_prop_time[4],
                                duration = v_duration[4]))
fun_vax <- approxfun(x = v_time_stop_vax, y = v_vax_rates, method = "linear")
fun_vax(0:100)

# Create list of parameters ----
l_parameters <- list(n_pop_size = n_pop_size,
                     n_inf   = n_inf, 
                     r_beta  = r_beta,
                     n_contacts_comm = n_contacts_comm,
                     r_sigma = r_sigma,
                     r_gamma = r_gamma,
                     n_hhsize  = n_hhsize,
                     n_contacts_hh = n_contacts_hh, 
                     r_tau   = r_tau, 
                     r_omega = r_omega,
                     r_vax_omega = r_vax_omega,
                     r_dx          = r_dx,
                     p_alpha_dx    = p_alpha_dx,
                     r_death       = r_death       ,
                     p_death_inf   = p_death_inf   ,
                     r_growth_rate = r_growth_rate ,
                     r_birth       = r_birth       ,
                     n_exp_states  = n_exp_states, 
                     n_inf_states  = n_inf_states,
                     times = times,
                     fun_npi = fun_npi,
                     fun_vax = fun_vax
)

## Verify NPI is working for a specific time ----
get_npi(n_time = 70, parameters = l_parameters)

# Run the model ----
sim_time <- system.time(
  l_out <- hh_mc_seir_out(parameters = l_parameters)
)
sim_time

# Plot epidemic outputs ----
show_MC_EI_model_results(l_out)
show_MC_SEIRV_model_results(l_out)

