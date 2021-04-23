rm(list = ls())

devtools::load_all()

library(tidyverse)
library(deSolve)
library(dplyr)
library(foreach)

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

# r_vax <- # TBD

r_dx          <- 0.1  # detection rate
p_alpha_dx    <- 0.20  # decrease in infectiousness
r_death       <- 0 #0.001 
p_death_inf   <- 0.02
r_growth_rate <- 1# 1.01 # population growth rate
r_birth       <- r_death*r_growth_rate

max_time <- 100

times <- seq(0, max_time, by = 1)

### NPIs
v_npi_times  <- c(0, 15, 16, 70, max_time+1)
# v_npi_levels <- c(1, 1, 0.2, 0.2, 0.2)
v_npi_levels <- c(1, 1, 0.01, 0.01, 0.01)
# v_npi_levels <- c(1, 1, 0.6, 0.6, 0.6)
# v_npi_levels <- c(1, 1, 1, 1, 1)

fun_npi <- approxfun(x = v_npi_times, y = v_npi_levels, method = "linear")

### Vax rates
### vaccination strategies
v_time_stop_vax <- c(0, 40, 41, max_time+1)
v_duration  <- diff(c(0, v_time_stop_vax))
# v_cum_prop_time <- c(0, 0.80, 0, 0)
v_cum_prop_time <- c(0, 0, 0, 0)

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
system.time(
  l_out <- hh_mc_seir_out(parameters = l_parameters)
)


# show_MC_EI_model_results(l_out)
show_MC_SEIRV_model_results(l_out)

v_hhsize  <- seq(1, 7)
v_num_exp <- seq(1, 3)
v_num_inf <- seq(1, 3)
v_r_beta  <- c(0.25, 0.35)
v_r_tau   <- c(0.40, 0.50)
v_r_omega <- c(0, 1/100, 1/50)
v_level_npi <- c(1, 0.8, 0.4)
v_vax_eff   <- c(1, 0.9, 0.5)
v_vax_prop  <- c(0, 0.3, 0.6, 0.9)
# n_hh_mod <- factorial(n_hhsize + (n_states-1))/(factorial(n_hhsize)*factorial(n_states-1)) #l_params_all

df_doe_mc_seirv <- expand.grid(n_hhsize = v_hhsize,
                               n_exp_states = v_num_exp,
                               n_inf_states = v_num_inf,
                               r_beta = v_r_beta,
                               r_tau = v_r_tau,
                               r_omega = v_r_omega,
                               level_npi = v_level_npi,
                               eff_vax = v_vax_eff,
                               vax_prop = v_vax_prop) %>%
  arrange(n_hhsize, n_exp_states, n_inf_states) %>%
  mutate(pid = row_number())
save(df_doe_mc_seirv, file = "data/df_doe_mc_seirv.RData")
saveRDS(df_doe_mc_seirv, file = "data/df_doe_mc_seirv.RDS")

#### Natural history DoE
df_doe_mc_seirv_nathist <- df_doe_mc_seirv %>% 
  filter(level_npi == 1,
         eff_vax   == 1,
         vax_prop  == 0)
save(df_doe_mc_seirv_nathist, file = "data/df_doe_mc_seirv_nathist.RData")
saveRDS(df_doe_mc_seirv_nathist, file = "data/df_doe_mc_seirv_nathist.RDS")

df_doe_mc_seirv_nathist_hhsize <- df_doe_mc_seirv_nathist %>% 
  filter(r_beta  == 0.25, 
         r_tau   == 0.40, 
         r_omega == 0)
save(df_doe_mc_seirv_nathist_hhsize, file = "data/df_doe_mc_seirv_nathist_hhsize.RData")
saveRDS(df_doe_mc_seirv_nathist_hhsize, file = "data/df_doe_mc_seirv_nathist_hhsize.RDS")

#### Control measures DoE
df_doe_mc_seirv_control <- df_doe_mc_seirv %>% 
  filter(n_hhsize %in% c(1, 3, 5), 
         # n_exp_states == 3, n_inf_states == 2,
         r_omega %in% c(0, 0.01, 0.02),
         r_tau %in% c(0.40, 0.50),
         r_beta %in% c(0.25, 0.35))
save(df_doe_mc_seirv_control, file = "data/df_doe_mc_seirv_control.RData")
saveRDS(df_doe_mc_seirv_control, file = "data/df_doe_mc_seirv_control.RDS")

