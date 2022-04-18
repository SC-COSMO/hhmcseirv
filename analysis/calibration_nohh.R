#### Analyze epidemic outputs ####

rm(list = ls())

# General setup ----
## Load packages ----
library(dplyr)
library(ggplot2)
library(ggpubr)
library(pBrackets)
library(patchwork)
library(stargazer)
library(grid)
library(png)

devtools::load_all()

library(tidyverse)
library(deSolve)
library(dplyr)
library(foreach)
# install.packages("remotes")
# remotes::install_github("coolbutuseless/ggpattern")
# library(ggpattern) # https://coolbutuseless.github.io/2020/04/01/introducing-ggpattern-pattern-fills-for-ggplot/
# source("R/02_decision_model_functions.R")

# Control measure outputs ----
## Load data ----
load(file = "output/df_output_doe_mc_seirv_all_control.RData")

## Wrangle data ----
# Rename variables for pretty plotting format
df_out_inf_all$Esize <- paste0("# of E compartments = ", df_out_inf_all$n_exp_states)
df_out_inf_all$Isize <- paste0("# of I compartments = ", df_out_inf_all$n_inf_states)
df_out_inf_all$`Household size` <- ordered(df_out_inf_all$n_hhsize, unique(df_out_inf_all$n_hhsize))
df_out_inf_all$`Household size labels` <- paste0("Household size = ", df_out_inf_all$n_hhsize)
# df_out_inf_all$n_hhsize <- ordered(df_out_inf_all$n_hhsize)
df_out_inf_all$`Vaccine effectiveness` <- scales::percent(df_out_inf_all$eff_vax)
df_out_inf_all$PropVax <- paste0("Proportion vaccinated = ", scales::percent(df_out_inf_all$vax_prop))
df_out_inf_all$EffVax  <- paste0("Vaccine effectiveness = ", scales::percent(df_out_inf_all$eff_vax))
df_out_inf_all$NPIeff <- paste0("NPI effectiveness = ", scales::percent(1-df_out_inf_all$level_npi))
df_out_inf_all$NPIeff_labels <- ordered(df_out_inf_all$NPIeff,
                                        unique(df_out_inf_all$NPIeff), c("No NPI", "NPI20", "NPI"))
df_out_inf_all$NPIeff_simple <- paste0(scales::percent(1-df_out_inf_all$level_npi))
df_out_inf_all$`Multicompartment structure` <- paste0("E=", 
                                                      df_out_inf_all$n_exp_states, 
                                                      ", I=", 
                                                      df_out_inf_all$n_inf_states)

### Subset to only no interventions ----
df_out_inf_noint <- df_out_inf_all %>%
  filter(vax_prop == 0,
         level_npi == 1)

df_out_inf_noint_hh <- df_out_inf_noint %>%
  filter(n_hhsize > 1, 
         time < 15)

df_out_inf_noint_nohh <- df_out_inf_noint_hh %>%
  filter(time == 0) %>%
  mutate(r_tau = 0,
          n_hhsize = 1)

run_hhmcseirv <- function(df_params, max_time = 14, nat_hist = TRUE){
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
  if(!nat_hist){ # If control measures
    ### NPIs
    v_npi_times  <- c(0, 9, 10, 70, max_time+1) 
    # v_npi_levels <- c(1, 1, 0.2, 0.2, 0.2)
    v_npi_levels <- c(1, 1, df_params$level_npi, df_params$level_npi, df_params$level_npi)
    # v_npi_levels <- c(1, 1, 0.6, 0.6, 0.6)
    # v_npi_levels <- c(1, 1, 1, 1, 1)
    
    fun_npi <- approxfun(x = v_npi_times, y = v_npi_levels, method = "linear")
    
    ### Vax rates
    ### vaccination strategies
    v_time_stop_vax <- c(0, 9, 10, 11, max_time+1)
    v_duration  <- diff(c(0, v_time_stop_vax))
    # v_cum_prop_time <- c(0, 0.80, 0, 0)
    v_cum_prop_time <- c(0, 0, df_params$vax_prop, 0, 0)
    
    v_vax_rates <- cbind(#daily_rate(cum_prop = v_cum_prop_time[1], # We need constant rate to get up to cumulative coverage
      #           duration = v_duration[1]),       # We need constant rate to get up to cumulative coverage
      0,
      daily_rate(cum_prop = v_cum_prop_time[2],
                 duration = v_duration[2]),
      daily_rate(cum_prop = v_cum_prop_time[3],
                 duration = v_duration[3]),
      daily_rate(cum_prop = v_cum_prop_time[4],
                 duration = v_duration[4]),
      0)
    fun_vax <- approxfun(x = v_time_stop_vax, y = v_vax_rates, method = "linear")
  } else{ # If Natural history
    ### NPIs
    # v_npi_times  <- c(0, 15, 16, 70, max_time+1)
    # v_npi_levels <- c(1, 1, 0.2, 0.2, 0.2)
    # v_npi_levels <- c(1, 1, df_params$level_npi, df_params$level_npi, df_params$level_npi)
    # v_npi_levels <- c(1, 1, 0.6, 0.6, 0.6)
    # v_npi_levels <- c(1, 1, 1, 1, 1)
    
    fun_npi <- approxfun(x = c(0:(max_time + 1)), y = rep(1, (max_time + 2)), method = "linear")
    
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
  }
  
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
  
  v_Inftot <- calc_inf_totals(l_out)$Inftot
  
  return(list(l_out = l_out,
              v_Inftot = v_Inftot))
}

gof <- function(r_beta, v_targets, df_params){
  df_params$r_beta <- r_beta
  out_model <- run_hhmcseirv(df_params)
  mse <- sum((v_targets - out_model$v_Inftot)^2)
  return(mse)
}

# 
# m_targets <- matrix(NA, 
#                     nrow = nrow(df_out_inf_noint_nohh), 
#                     ncol = 15)
# l_out_optim <- vector(mode = "list", length = nrow(df_out_inf_noint_nohh))
# 
# i <- 1
# for (n_pid in unique(df_out_inf_noint_nohh$pid)[1]) {
#   m_targets[i, ] <- df_out_inf_noint_hh %>%
#     filter(pid == n_pid) %>%
#     select(Inftot) %>%
#     t()
#   v_targets <- m_targets[i, ]
#   df_params <- df_out_inf_noint_nohh %>%
#     filter(pid == n_pid)
#   
#   l_out_optim[[i]] <- optim(par = c(r_beta = df_params$r_beta), 
#                             fn = gof, lower = 0, upper = 0.8, 
#                             method = "Brent",
#                             v_targets = v_targets, df_params = df_params)
#   print(i)
#   i <- i + 1
# }
# gof(r_beta = 0.3, v_targets = v_targets, df_params = df_params)

no_cores <- parallel::detectCores() - 1
cl <- parallel::makeForkCluster(no_cores) 
doParallel::registerDoParallel(cl)

l_out_optim <- foreach::foreach(n_pid = unique(df_out_inf_noint_nohh$pid))  %dopar% {
  v_targets <- df_out_inf_noint_hh %>%
    filter(pid == n_pid) %>%
    select(Inftot) %>%
    t()
  df_params <- df_out_inf_noint_nohh %>%
    filter(pid == n_pid)
  
  l_out_optim_single <- optim(par = c(r_beta = df_params$r_beta), 
                              fn = gof, lower = 0, upper = 0.8, 
                              method = "Brent",
                              v_targets = v_targets, df_params = df_params)
  save(l_out_optim_single, file = paste("output/calibration/f_", n_pid, ".RData"))
  l_out_optim_single
}
parallel::stopCluster(cl)

l_out_optim <- vector(mode = "list", length = nrow(df_out_inf_noint_nohh))
df_out_optim <- c()
i <- 1
for(n_pid in unique(df_out_inf_noint_nohh$pid)){
  load(paste("output/calibration/f_", n_pid, ".RData")) 
  l_out_optim[[i]] <- l_out_optim_single
  i <- i + 1
  df_out_optim_single <- data.frame(pid = n_pid,
                                     r_beta = l_out_optim_single$par,
                                     MSE = l_out_optim_single$value,
                                     counts = l_out_optim_single$counts[1],
                                     convergence = l_out_optim_single$convergence)
  df_out_optim <- rbind.data.frame(df_out_optim,
                                   df_out_optim_single)
}
save(df_out_optim, file = "output/df_out_optim.RData")
# hist(df_out_optim$r_beta)

load(file = "output/df_out_optim.RData")

df_out_inf_hh_gt1 <- df_out_inf_all %>%
  filter(n_hhsize > 1, 
         time == 0)

no_cores <- parallel::detectCores() - 2
cl <- parallel::makeForkCluster(no_cores) 
doParallel::registerDoParallel(cl)

l_out_projection <- foreach::foreach(n_pid = unique(df_out_optim$pid))  %dopar% {
  temp_r_beta <- df_out_optim %>%
    filter(pid == n_pid) %>%
    select(r_beta) %>%
    as.numeric()
  
  df_params_temp <- df_out_inf_hh_gt1 %>%
    filter(pid == n_pid)
  df_temp <- df_out_inf_hh_gt1 %>%
    filter(n_hhsize == df_params_temp$n_hhsize,
           n_exp_states == df_params_temp$n_exp_states,
           n_inf_states == df_params_temp$n_inf_states,
           r_beta == df_params_temp$r_beta,
           r_tau == df_params_temp$r_tau,
           r_omega ==df_params_temp$r_omega)  %>%
    select(n_hhsize:pid)
  df_temp$n_hhsize <- 1
  df_temp$r_tau <- 0
  df_temp$r_beta <- temp_r_beta
  df_out_projection_single <- c()
  for(i in 1:nrow(df_temp)){
    df_params <- df_temp[i, ]
    l_out_model <- run_hhmcseirv(df_params = df_params, 
                                 max_time = 100, 
                                 nat_hist = FALSE)
    df_out_projection_single <- bind_rows(df_out_projection_single,
                                          data.frame(df_params, 
                                                     l_out_model$l_out$df_out_hh_mc_seir))
  }
  df_out_projection_single
}
parallel::stopCluster(cl)
save(l_out_projection, 
     file = "output/projection_hh_gt1/l_out_projection.RData")
