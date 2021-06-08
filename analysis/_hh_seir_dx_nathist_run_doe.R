rm(list = ls())

devtools::load_all()

# install.packages("diffeqr")
# library(diffeqr)
library(tidyverse)
library(deSolve)
library(dplyr)
library(foreach)

# Read the location as an argument
args <- commandArgs(trailingOnly = TRUE)
n_pid <- ifelse(!is.na(args[1]), args[1], "1") 
n_pid <- as.numeric(n_pid)
print(n_pid)

df_doe_mc_seirv <- readRDS("data/df_doe_mc_seirv.RDS")

df_params <- df_doe_mc_seirv %>%
  filter(pid == n_pid)


n_pop_size <- 1000000
n_inf   <- 1
n_contacts_comm <- 10 
n_contacts_hh <- 4
r_vax_omega <- 0 
r_sigma <- 0.35
r_gamma <- 0.2
r_dx          <- 0.1  # detection rate
p_alpha_dx    <- 0.20  # decrease in infectiousness
r_death       <- 0 #0.001 
p_death_inf   <- 0.02
r_growth_rate <- 1# 1.01 # population growth rate
r_birth       <- r_death*r_growth_rate

max_time <- 200

times <- seq(0, max_time, by = 1)

### DoE parameters
n_hhsize <- df_params$n_hhsize
r_beta   <- df_params$r_beta
r_tau    <- df_params$r_tau
r_omega  <- df_params$r_omega
eff_vax  <- df_params$eff_vax

n_exp_states <- df_params$n_exp_states
n_inf_states <- df_params$n_inf_states

# r_vax <- # TBD

### NPIs
v_npi_times  <- c(0, 15, 16, 70, max_time+1)
# v_npi_levels <- c(1, 1, 0.2, 0.2, 0.2)
v_npi_levels <- c(1, 1, df_params$level_npi, df_params$level_npi, df_params$level_npi)
# v_npi_levels <- c(1, 1, 0.6, 0.6, 0.6)
# v_npi_levels <- c(1, 1, 1, 1, 1)

fun_npi <- approxfun(x = v_npi_times, y = v_npi_levels, method = "linear")

### Vax rates
### vaccination strategies
v_time_stop_vax <- c(0, 40, 41, max_time+1)
v_duration  <- diff(c(0, v_time_stop_vax))
# v_cum_prop_time <- c(0, 0.80, 0, 0)
v_cum_prop_time <- c(0, df_params$vax_prop, 0, 0)

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
# list2env(l_parameters, envir = .GlobalEnv)

get_npi(n_time = 70, parameters = l_parameters)

## Test run
sim_time <- system.time(
  l_out <- hh_mc_seir_out(parameters = l_parameters)
)
l_out$l_params_all$sim_time <- sim_time[3]

save(l_out, file = paste0("output/output_doe_mc_seirv_", n_pid,".RData"))
