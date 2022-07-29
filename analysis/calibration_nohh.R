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
library(xtable)
library(grid)
library(png)
library(multcomp)

devtools::load_all()

library(tidyverse)
library(deSolve)
library(dplyr)
library(foreach)
# install.packages("remotes")
# remotes::install_github("coolbutuseless/ggpattern")
# library(ggpattern) # https://coolbutuseless.github.io/2020/04/01/introducing-ggpattern-pattern-fills-for-ggplot/
# source("R/02_decision_model_functions.R")

# Control measure outputs with Household structure----
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
df_out_inf_all$NPIeff <- paste0("NPI effectiveness = ", scales::percent(1 - df_out_inf_all$level_npi))
df_out_inf_all$NPIeff_labels <- ordered(df_out_inf_all$NPIeff,
                                        unique(df_out_inf_all$NPIeff), c("No NPI", "NPI20", "NPI60"))
df_out_inf_all$NPIeff_simple <- paste0(scales::percent(1 - df_out_inf_all$level_npi))
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
                       fun_vax = fun_vax,
                       eff_vax = eff_vax
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


# Generate projections with calibrated betas in the absence of HH structure ----
load(file = "output/df_out_optim.RData")

  
df_out_inf_hh_gt1 <- df_out_inf_all %>%
  filter(n_hhsize > 1, 
         time == 0)

no_cores <- parallel::detectCores() - 2
cl <- parallel::makeForkCluster(no_cores) 
doParallel::registerDoParallel(cl)

l_out_projection <- foreach::foreach(n_pid = unique(df_out_optim$pid))  %dopar% {
  temp_r_beta <- df_out_optim %>% # n_pid <- unique(df_out_optim$pid)[1]
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
  for(i in 1:nrow(df_temp)){ # i <- 1
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

# Compute bias Delta(r_beta, r_beta_nhh) ----
load(file = "output/projection_hh_gt1/l_out_projection.RData")

df_out_inf_hh_gt1_noint <- df_out_inf_hh_gt1 %>% 
  filter(vax_prop == 0,
         level_npi == 1)

df_out_optim_bias <- df_out_optim %>%
  rename(r_beta_nhh = r_beta) %>%
  right_join(df_out_inf_hh_gt1_noint, by = "pid") %>%
  mutate(delta_r_beta = r_beta_nhh - r_beta)

ggplot(df_out_optim_bias, aes(x = r_omega, y = r_beta_nhh - r_beta, color = as.factor(n_hhsize))) +
  facet_grid(Esize ~ Isize) +
  geom_point()

## Metaregression on bias Delta(r_beta, r_beta_nhh) ----
lm_out_delta_r_beta <- lm(delta_r_beta ~ n_exp_states*n_inf_states + 
     n_exp_states*n_hhsize + 
     n_inf_states*n_hhsize +
     r_tau + r_beta, data = df_out_optim_bias)
summary(lm_out_delta_r_beta)

# Compute summary projections from No HH projections ----
n_proj_nhh <- length(l_out_projection)
df_out_inf_all_nhh <- c()
for (i in 1:n_proj_nhh) { # i <- 2
  temp <- l_out_projection[[i]]
  v_names_exp    <- paste("E", letters[seq(1, temp$n_exp_states[1])], sep = "")
  v_names_inf    <- paste("I", letters[seq(1, temp$n_inf_states[1])], sep = "")
  v_names_inf_dx <- paste("IDX", letters[seq(1, temp$n_inf_states[1])], sep = "")
  
  df_temp <- data.frame(pid = temp$pid, 
                        level_npi = temp$level_npi,
                        eff_vax = temp$eff_vax,
                        vax_prop = temp$vax_prop,
                        time = temp$time,
                        n_hhsize = temp$n_hhsize,
                        r_beta = temp$r_beta,
                        r_tau = temp$r_tau,
                        r_omega = temp$r_omega,
                        n_exp_states = temp$n_exp_states,
                        n_inf_states = temp$n_inf_states,
                        # Exptot_nhh = rowSums(temp[, v_names_exp,drop=FALSE]),
                        InfNoDX_nhh = rowSums(temp[, v_names_inf,drop=FALSE]),
                        Inftot_nhh = rowSums(temp[, c(v_names_inf,
                                                    v_names_inf_dx), drop=FALSE]),
                        check.names = FALSE)
  df_out_inf_all_nhh <- bind_rows(df_out_inf_all_nhh,
                                  df_temp)
}
save(df_out_inf_all_nhh, 
     file = "output/projection_hh_gt1/df_out_inf_all_nhh.RData")

load(file = "output/projection_hh_gt1/df_out_inf_all_nhh.RData")

# Rename variables for pretty plotting format
df_out_inf_all_nhh$Esize <- paste0("# of E compartments = ", df_out_inf_all_nhh$n_exp_states)
df_out_inf_all_nhh$Isize <- paste0("# of I compartments = ", df_out_inf_all_nhh$n_inf_states)
df_out_inf_all_nhh$`Household size` <- ordered(df_out_inf_all_nhh$n_hhsize, unique(df_out_inf_all_nhh$n_hhsize))
df_out_inf_all_nhh$`Household size labels` <- paste0("Household size = ", df_out_inf_all_nhh$n_hhsize)
# df_out_inf_all_nhh$n_hhsize <- ordered(df_out_inf_all_nhh$n_hhsize)
df_out_inf_all_nhh$`Vaccine effectiveness` <- scales::percent(df_out_inf_all_nhh$eff_vax)
df_out_inf_all_nhh$PropVax <- paste0("Proportion vaccinated = ", scales::percent(df_out_inf_all_nhh$vax_prop))
df_out_inf_all_nhh$EffVax  <- paste0("Vaccine effectiveness = ", scales::percent(df_out_inf_all_nhh$eff_vax))
df_out_inf_all_nhh$NPIeff <- paste0("NPI effectiveness = ", scales::percent(1 - df_out_inf_all_nhh$level_npi))
df_out_inf_all_nhh$NPIeff_labels <- ordered(df_out_inf_all_nhh$NPIeff,
                                                  unique(df_out_inf_all_nhh$NPIeff), c("No NPI", "NPI20", "NPI60"))
df_out_inf_all_nhh$NPIeff_simple <- paste0(scales::percent(1 - df_out_inf_all_nhh$level_npi))
df_out_inf_all_nhh$`Multicompartment structure` <- paste0("E=", 
                                                          df_out_inf_all_nhh$n_exp_states, 
                                                          ", I=", 
                                                          df_out_inf_all_nhh$n_inf_states)

df_out_inf_all_nhh_noint <- df_out_inf_all_nhh %>% 
  filter(vax_prop == 0,
         level_npi == 1)

df_out_inf_all_nhh_noint_summ <- df_out_inf_all_nhh_noint %>%
  group_by(pid) %>%
  mutate(# Find the t at which I(t) is at its max
    max_Inftot_nhh = max(Inftot_nhh),
    max_Inftot_nhh_time = time[which.max(Inftot_nhh)],
    # Find times at which I(t) is at a percentage of its max
    p05_Inftot_nhh_time = max(time[which(Inftot_nhh <= max_Inftot_nhh*0.05 & time < max_Inftot_nhh_time)]),
    p10_Inftot_nhh_time = max(time[which(Inftot_nhh <= max_Inftot_nhh*0.10 & time < max_Inftot_nhh_time)]),
    p25_Inftot_nhh_time = max(time[which(Inftot_nhh <= max_Inftot_nhh*0.25 & time < max_Inftot_nhh_time)]),
    p50_Inftot_nhh_time = max(time[which(Inftot_nhh <= max_Inftot_nhh*0.50 & time < max_Inftot_nhh_time)]),
    # Find times at which x number of IDX(t) are seen
    IDX500_nhh_time = max(time[which((Inftot_nhh - InfNoDX_nhh) <= 500 & time < max_Inftot_nhh_time)]),
    IDX100_nhh_time = max(time[which((Inftot_nhh - InfNoDX_nhh) <= 100 & time < max_Inftot_nhh_time)]),
    CumInfTot_nhh = sum(Inftot_nhh)) %>%
  slice_head() %>%
  ungroup()

# Control measures outputs with Household structure ----
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
df_out_inf_all$NPIeff <- paste0("NPI effectiveness = ", scales::percent(1 - df_out_inf_all$level_npi))
df_out_inf_all$NPIeff_labels <- ordered(df_out_inf_all$NPIeff,
                                        unique(df_out_inf_all$NPIeff), c("No NPI", "NPI20", "NPI60"))
df_out_inf_all$NPIeff_simple <- paste0(scales::percent(1 - df_out_inf_all$level_npi))
df_out_inf_all$`Multicompartment structure` <- paste0("E=", 
                                                      df_out_inf_all$n_exp_states, 
                                                      ", I=", 
                                                      df_out_inf_all$n_inf_states)

## Summarize output ----
df_out_inf_all_summ <- df_out_inf_all %>%
  group_by(pid) %>%
  mutate(# Find the t at which I(t) is at its max
    ref_category = ifelse(n_hhsize == 1 & n_exp_states == 1 & n_inf_states ==1, 
                          1, 0), # Define HH=1, E=1, I=1 as reference category
    max_Inftot = max(Inftot),
    max_Inftot_time = time[which.max(Inftot)],
    # Find times at which I(t) is at a percentage of its max
    p05_Inftot_time = max(time[which(Inftot <= max_Inftot*0.05 & time < max_Inftot_time)]),
    p10_Inftot_time = max(time[which(Inftot <= max_Inftot*0.10 & time < max_Inftot_time)]),
    p25_Inftot_time = max(time[which(Inftot <= max_Inftot*0.25 & time < max_Inftot_time)]),
    p50_Inftot_time = max(time[which(Inftot <= max_Inftot*0.50 & time < max_Inftot_time)]),
    # Find times at which x number of IDX(t) are seen
    IDX500_time = max(time[which((Inftot-InfNoDX) <= 500 & time < max_Inftot_time)]),
    IDX100_time = max(time[which((Inftot-InfNoDX) <= 100 & time < max_Inftot_time)]),
    CumInfTot = sum(Inftot)) %>%
  slice_head() %>%
  ungroup()

# Combine no HH structure projections with natural history projections ----
df_out_inf_all_biased_noint_summ <- df_out_inf_all_nhh_noint_summ %>%
  dplyr::select(-n_hhsize, -r_beta, -r_tau, -`Household size`, -`Household size labels`) %>%
  right_join(df_out_inf_all_summ %>% 
               filter(n_hhsize > 1, 
                      vax_prop == 0,
                      level_npi == 1), 
             by = "pid") %>%
  mutate(max_Inftot_diff = max_Inftot_nhh - max_Inftot,
         max_Inftot_time_diff = max_Inftot_nhh_time - max_Inftot_time,
         CumInfTot_diff = CumInfTot_nhh - CumInfTot,
         # Percentage changes
         max_Inftot_diff_perc = (max_Inftot_diff/max_Inftot)*100,
         max_Inftot_time_diff_perc = (max_Inftot_time_diff/max_Inftot_time)*100,
         CumInfTot_diff_perc = (CumInfTot_diff/CumInfTot)*100)
df_out_inf_all_biased_noint_summ$HH_5 <- df_out_inf_all_biased_noint_summ$n_hhsize == 5

# Plot epidemic curves HH vs no HH structure ----
v_fig_pid <- df_out_inf_all_summ %>% filter(n_hhsize == 3, 
                                            eff_vax == 1 & 
                               vax_prop == 0,
                               level_npi == 1,
                               r_beta == 0.25 & r_tau == 0.50 & 
                               r_omega == 0.000 & time <= 60 &
                               `Multicompartment structure` %in% c("E=1, I=1", 
                                                                   "E=3, I=3", 
                                                                   "E=3, I=1", 
                                                                   "E=1, I=3") 
                               ) %>%
  dplyr::select(pid)

df_fig_nohh_vs_hh <-  bind_rows(df_out_inf_all_nhh_noint %>% 
                                  mutate(Structure = "No HH") %>%
                                  mutate(`Household size` = ordered(1, 
                                                                    levels = c(1, 3, 5)),
                                         `NPIeff_labels` = ordered("No NPI", 
                                                                   levels = c("No NPI", "NPI20", "NPI60"))) %>%
                                  rename(#Exptot = Exptot_nhh ,
                                         InfNoDX = InfNoDX_nhh,
                                         Inftot = Inftot_nhh), 
                                df_out_inf_all %>% 
                                  mutate(Structure = "HH")) %>% 
  filter(pid %in% as.matrix(v_fig_pid))



gg_epidemic_curve_nohh_vs_hh <- ggplot(df_fig_nohh_vs_hh, 
                                                aes(x = time, y = Inftot/10e6, 
                                                    color = Structure)) + # 
  geom_line(size = 1.1) +
  # geom_segment(aes(x = max_Inftot_time_E1_I1_hh1$max_Inftot_time, 
  #                  xend = max_Inftot_time_E1_I1_hh1$max_Inftot_time, 
  #                  y = 0, 
  #                  yend = max_Inftot_time_E1_I1_hh1$max_Inftot),
  #              linetype = "dashed", color = "blue") +
  # geom_segment(aes(x = max_Inftot_time_E1_I1_hh3$max_Inftot_time, 
  #                  xend = max_Inftot_time_E1_I1_hh3$max_Inftot_time, 
  #                  y = 0, 
  #                  yend = max_Inftot_time_E1_I1_hh3$max_Inftot),
  #              linetype = "dashed", color = "red") +
  # geom_point(aes(x = max_Inftot_time_E1_I1_hh1$max_Inftot_time, 
  #                y = max_Inftot_time_E1_I1_hh1$max_Inftot),
  #            color = "blue", size = 2.5) +
  # geom_point(aes(x = max_Inftot_time_E1_I1_hh3$max_Inftot_time, 
  #                y = max_Inftot_time_E1_I1_hh3$max_Inftot),
  #            color = "red", size = 2.5) +
  # geom_segment(aes(x = x0, xend = max_Inftot_time,
  #                  y = y0, yend = max_Inftot,
  #                  color = `Household size`),
  #              linejoin = "mitre",
  #              linetype = "dashed", arrow = arrow(length = unit(0.15, "inches"),
  #                                                 type = "closed"), 
  #              show.legend = FALSE) +
  # geom_vline(xintercept = as.numeric(max_Inftot_time_E1_I1_hh1), linetype = "dashed", color = "blue") +
  # geom_vline(xintercept = as.numeric(max_Inftot_time_E1_I1_hh3), linetype = "dashed", color = "red") +
  annotate("text",
           label  = "Calibration",
           x      = 8,
           y      = 0.062,
           size   = 6,
           colour = "#393F3E") +
  annotate("text",
           label  = "Projection",
           x      = 22,
           y      = 0.062,
           size   = 6,
           colour = "#393F3E") +
  geom_vline(xintercept = 15, linetype = "dashed", color = "black") +
  facet_grid(Esize ~ Isize) +
  scale_x_continuous(limits = c(0, 50)) +
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1.0), limits = c(0, 0.065)) +
  # scale_color_grey(start = 0.2, end = 0.6) +
  scale_color_manual(values = c("No HH" = "blue", "HH" = "red")) +
  xlab("Time") +
  ylab("Infected population (% of total population)") +
  guides(color = guide_legend(nrow = 2), 
         linetype = guide_legend(nrow = 1)) +
  theme_bw(base_size = 20) +
  theme(strip.background = element_rect(colour = "white", fill = "white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        # legend.position = c(""),
        legend.position = c(0.85, 0.90),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key = element_blank())
gg_epidemic_curve_nohh_vs_hh
ggsave(plot = gg_epidemic_curve_nohh_vs_hh, 
       filename = "figs/Paper/Fig2_natural_history_curves_NoHH_HH.pdf", 
       width = 12, height = 8)

# Metaregression ----
## Without interactions
fit_hh_peak_time_bias_rel_nhh <- lm(max_Inftot_time_diff_perc ~ n_exp_states.x + 
                                      n_inf_states.x + 
                                      HH_5 + 
                                      r_omega.x +
                                      r_tau + r_beta, 
                                    data = df_out_inf_all_biased_noint_summ %>% 
                                      filter(max_Inftot_time_diff_perc != Inf))
summary(fit_hh_peak_time_bias_rel_nhh)

fit_hh_peak_size_bias_rel_nhh <- lm(max_Inftot_diff_perc ~ n_exp_states.x + 
                                      n_inf_states.x + 
                                      HH_5 + 
                                      r_omega.x +
                                      r_tau + r_beta,
                                    data = df_out_inf_all_biased_noint_summ)
summary(fit_hh_peak_size_bias_rel_nhh)

fit_hh_epidemic_size_rel_nhh <- lm(CumInfTot_diff_perc ~ n_exp_states.x + 
                                     n_inf_states.x + 
                                     HH_5 + 
                                     r_omega.x +
                                     r_tau + r_beta, 
                                   data = df_out_inf_all_biased_noint_summ)
summary(fit_hh_epidemic_size_rel_nhh)

# Linear combination coefficient tests
l_lincomm_rel_nhh <- vector(mode = "list", length = 3)
l_lincomm_rel_nhh[[1]] <- summary(multcomp::glht(fit_hh_peak_time_bias_rel_nhh, linfct = c("(Intercept) + HH_5TRUE = 0")))
l_lincomm_rel_nhh[[2]] <- summary(multcomp::glht(fit_hh_peak_size_bias_rel_nhh, linfct = c("(Intercept) + HH_5TRUE = 0")))
l_lincomm_rel_nhh[[3]] <- summary(multcomp::glht(fit_hh_epidemic_size_rel_nhh, linfct = c("(Intercept) + HH_5TRUE = 0")))

df_lincomm_rel_nhh <- data.frame(`Epidemic outcome` = c("Peak time", "Peak size", "Epidemic size"),
                                 coef = as.numeric(c(l_lincomm_rel_nhh[[1]]$test$coefficients, 
                                                     l_lincomm_rel_nhh[[2]]$test$coefficients,
                                                     l_lincomm_rel_nhh[[3]]$test$coefficients)), 
                                 pval = c(l_lincomm_rel_nhh[[1]]$test$pvalues[1], 
                                          l_lincomm_rel_nhh[[2]]$test$pvalues[1],
                                          l_lincomm_rel_nhh[[3]]$test$pvalues[1]), 
                                 check.names = F)
df_lincomm_rel_nhh$Interaction <- "No"
df_lincomm_rel_nhh
xtable::xtable(df_lincomm_rel_nhh, )

## With interactions
fit_hh_peak_time_bias_rel_nhh_int <- lm(max_Inftot_time_diff_perc ~ #n_exp_states.x*n_inf_states.x + 
                                      n_exp_states.x*HH_5 + 
                                      n_inf_states.x*HH_5 + 
                                        r_omega.x +
                                      r_tau + r_beta, 
                                    data = df_out_inf_all_biased_noint_summ %>% 
                                      filter(max_Inftot_time_diff_perc != Inf))
summary(fit_hh_peak_time_bias_rel_nhh_int)

fit_hh_peak_size_bias_rel_nhh_int <- lm(max_Inftot_diff_perc ~ # n_exp_states.x*n_inf_states.x + 
                                          n_exp_states.x*HH_5 + 
                                          n_inf_states.x*HH_5 + 
                                          r_omega.x +
                                      r_tau + r_beta,
                                    data = df_out_inf_all_biased_noint_summ)
summary(fit_hh_peak_size_bias_rel_nhh_int)

fit_hh_epidemic_size_rel_nhh_int <- lm(CumInfTot_diff_perc ~ #n_exp_states.x*n_inf_states.x + 
                                         n_exp_states.x*HH_5 + 
                                         n_inf_states.x*HH_5 + 
                                         r_omega.x +
                                     r_tau + r_beta, 
                                   data = df_out_inf_all_biased_noint_summ)
summary(fit_hh_epidemic_size_rel_nhh_int)

# Linear combination coefficient tests
l_lincomm_rel_nhh_int <- vector(mode = "list", length = 3)
l_lincomm_rel_nhh_int[[1]] <- summary(multcomp::glht(fit_hh_peak_time_bias_rel_nhh_int, linfct = c("(Intercept) + HH_5TRUE + n_exp_states.x:HH_5TRUE + HH_5TRUE:n_inf_states.x = 0")))
l_lincomm_rel_nhh_int[[2]] <- summary(multcomp::glht(fit_hh_peak_size_bias_rel_nhh_int, linfct = c("(Intercept) + HH_5TRUE + n_exp_states.x:HH_5TRUE + HH_5TRUE:n_inf_states.x = 0")))
l_lincomm_rel_nhh_int[[3]] <- summary(multcomp::glht(fit_hh_epidemic_size_rel_nhh_int, linfct = c("(Intercept) + HH_5TRUE + n_exp_states.x:HH_5TRUE + HH_5TRUE:n_inf_states.x = 0")))

df_lincomm_rel_nhh_int <- data.frame(`Epidemic outcome` = c("Peak time", "Peak size", "Epidemic size"),
                                 coef = as.numeric(c(l_lincomm_rel_nhh_int[[1]]$test$coefficients, 
                                                     l_lincomm_rel_nhh_int[[2]]$test$coefficients,
                                                     l_lincomm_rel_nhh_int[[3]]$test$coefficients)), 
                                 pval = c(l_lincomm_rel_nhh_int[[1]]$test$pvalues[1], 
                                          l_lincomm_rel_nhh_int[[2]]$test$pvalues[1],
                                          l_lincomm_rel_nhh_int[[3]]$test$pvalues[1]), 
                                 check.names = F)
df_lincomm_rel_nhh_int$Interaction <- "Yes"

xtable::xtable(df_lincomm_rel_nhh_int)
# print(xtable(df_lincomm_rel_nhh_comb, 
#              caption = "Linear combinations estimates for HH5."),
#       include.rownames = FALSE)

df_lincomm_rel_nhh_comb <- bind_rows(df_lincomm_rel_nhh, df_lincomm_rel_nhh_int) %>%
  group_by(`Epidemic outcome`) %>%
  arrange(desc(`Epidemic outcome`), Interaction)
df_lincomm_rel_nhh_comb

v_pval_rel_nhh <- df_lincomm_rel_nhh_comb$pval
v_star_rel_nhh <- cut(v_pval_rel_nhh,
              breaks = c(-Inf, 0.01, 0.05, 0.1, Inf), 
              labels = c("$^{***}$", "$^{**}$", "$^{*}$", ""))

v_coef_rel_nhh <- as.character(round(df_lincomm_rel_nhh_comb$coef, 3))
v_coef_sg_rel_nhh <- paste0(v_coef_rel_nhh, v_star_rel_nhh)

v_stargazer_notes <- c("HH3 describes how  much the outcome differs if household structure is excluded",
                       "from an otherwise similar epidemic where the true household size is 3.",
                       "HH5 describes the incremental difference in the outcome if household structure is excluded",
                       "from an otherwise similar epidemic where the true household size is 5 instead of 3.",
                       "For the models, the total effect of exclusion of household structure involves interactions with other terms.",
                       "For example, E*HH5 describes how the incremental difference depends upon the number of exposed",
                       "compartments when the true household size is 5 instead of 3.",
                       "We estimated the magnitude of the linear combination of the relevant coefficients and tested their",
                       "significance.")
stargazer(fit_hh_peak_time_bias_rel_nhh, fit_hh_peak_time_bias_rel_nhh_int,
          fit_hh_peak_size_bias_rel_nhh, fit_hh_peak_size_bias_rel_nhh_int, 
          fit_hh_epidemic_size_rel_nhh, fit_hh_epidemic_size_rel_nhh_int, 
          type = "text", # latex 
          title = "Metaregression estimates on percentage differences in outcomes of excluding household structure",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          label = "tab:metaregepi",
          order = c(9, 3, 1, 2, 5, 6, 4, 7, 8),
          covariate.labels = c("Household size 3 (HH3)", "Household size 5 (HH5)",
                               "E", "I", "$\\tau$", "$\\beta$", "$\\omega$",
                               # "E*I",
                               "E*HH5", "I*HH5"),
          add.lines = list(c("HH3 + HH5", 
                             v_coef_sg_rel_nhh)),
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"), 
          notes = v_stargazer_notes, 
          notes.append = TRUE,
          notes.align = "l",
          align = TRUE) 
stargazer(fit_hh_peak_time_bias_rel_nhh, fit_hh_peak_time_bias_rel_nhh_int,
          fit_hh_peak_size_bias_rel_nhh, fit_hh_peak_size_bias_rel_nhh_int, 
          fit_hh_epidemic_size_rel_nhh, fit_hh_epidemic_size_rel_nhh_int,
          type = "latex",  
          out =  "output/Table3_metaregression_bias.tex",
          title = "Metaregression estimates on percentage differences in outcomes of excluding household structure",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          label = "tab:metaregepi",
          order = c(9, 3, 1, 2, 5, 6, 4, 7, 8),
          covariate.labels = c("Household size 3 (HH3)", "Household size 5 (HH5)",
                               "E", "I", "$\\tau$", "$\\beta$", "$\\omega$",
                               # "E*I",
                               "E*HH5", "I*HH5"),
          add.lines = list(c("HH3 + HH5", 
                             v_coef_sg_rel_nhh)),
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          notes = v_stargazer_notes, 
          notes.append = TRUE,
          notes.align = "l",
          align = TRUE) 

# Control measure outputs on No HH dataset----
df_out_inf_all_nhh

## Summarize output ----
df_out_inf_all_control_nhh_summ <- df_out_inf_all_nhh %>%
  group_by(pid) %>% # ,  n_exp_states, n_inf_states
  mutate(# Find the t at which I(t) is at its max
    # ref_category = ifelse(n_hhsize == 1, 
    #                       1, 0), # Define HH=1 as reference category
    # ref_category_bias = ifelse(n_hhsize == 1 & n_exp_states == 1 & n_inf_states == 1.0, 
    #                            1, 0), # n_hhsize = 1, E = 1, I =1 as reference category to estimate bias of control measures' effect
    ref_category_nathist = ifelse(level_npi == 1 & vax_prop == 0 & eff_vax == 1.0, 
                                  1, 0), # Natural history as reference category to estimate control measures' effectiveness on other parameters
    max_Inftot_nhh = max(Inftot_nhh),
    max_Inftot_nhh_time = time[which.max(Inftot_nhh)],
    # Find times at which I(t) is at a percentage of its max
    p05_Inftot_nhh_time = max(time[which(Inftot_nhh <= max_Inftot_nhh*0.05 & time < max_Inftot_nhh_time)]),
    p10_Inftot_nhh_time = max(time[which(Inftot_nhh <= max_Inftot_nhh*0.10 & time < max_Inftot_nhh_time)]),
    p25_Inftot_nhh_time = max(time[which(Inftot_nhh <= max_Inftot_nhh*0.25 & time < max_Inftot_nhh_time)]),
    p50_Inftot_nhh_time = max(time[which(Inftot_nhh <= max_Inftot_nhh*0.50 & time < max_Inftot_nhh_time)]),
    # Find times at which x number of IDX(t) are seen
    IDX500_nhh_time = max(time[which((Inftot_nhh - InfNoDX_nhh) <= 500 & time < max_Inftot_nhh_time)]),
    IDX100_nhh_time = max(time[which((Inftot_nhh - InfNoDX_nhh) <= 100 & time < max_Inftot_nhh_time)]),
    CumInfTot_nhh = sum(Inftot_nhh)) %>%
  # arrange(-ref_category_nathist) %>%
  slice_head() %>%
  ungroup() %>%
  group_by(n_hhsize, r_beta, r_tau, r_omega, n_exp_states, n_inf_states) %>%
  mutate(
    pid_nathist = max(pid*ref_category_nathist)
    ) %>%
  ungroup()

## Calculate  differential outcomes compared to reference category (Group by parameter set) ----
# Create a data.frame with only reference category (ie, hhsize = 1)
df_out_inf_all_nathist_nhh_summ <- df_out_inf_all_control_nhh_summ %>%
  filter(ref_category_nathist == 1) %>%
  rename(max_Inftot_nhh_nathist = max_Inftot_nhh,
         max_Inftot_nhh_time_nathist = max_Inftot_nhh_time,
         CumInfTot_nhh_nathist = CumInfTot_nhh) %>%
  dplyr::select(pid, pid_nathist, n_hhsize, n_exp_states, n_inf_states,
         r_beta, r_tau, r_omega, 
         max_Inftot_nhh_nathist, max_Inftot_nhh_time_nathist, CumInfTot_nhh_nathist)

# Compute differential outcomes compared to reference category (Group by parameter set)
df_out_inf_all_control_nathist_nhh_summ <- df_out_inf_all_control_nhh_summ %>%
  left_join(df_out_inf_all_nathist_nhh_summ, by = c("pid_nathist")) %>%
  mutate(max_Inftot_diff_nhh = max_Inftot_nhh - max_Inftot_nhh_nathist,
         max_Inftot_time_diff_nhh = max_Inftot_nhh_time - max_Inftot_nhh_time_nathist,
         CumInfTot_diff_nhh = CumInfTot_nhh - CumInfTot_nhh_nathist,
         # Percentage changes
         max_Inftot_diff_perc_nhh = (max_Inftot_diff_nhh/max_Inftot_nhh_nathist)*100,
         max_Inftot_time_diff_perc_nhh = (max_Inftot_time_diff_nhh/max_Inftot_nhh_time_nathist)*100,
         CumInfTot_diff_perc_nhh = (CumInfTot_diff_nhh/CumInfTot_nhh_nathist)*100) %>%
  rename(pid = pid.x)


## Wrangle data for dataset with Household structure
### Summarize output ----
df_out_inf_all_control_summ <- df_out_inf_all %>%
  filter(n_hhsize > 1) %>%
  group_by(pid) %>% # ,  n_exp_states, n_inf_states
  mutate(# Find the t at which I(t) is at its max
    ref_category_nathist = ifelse(level_npi == 1 & vax_prop == 0 & eff_vax == 1.0, 
                                  1, 0), # Natural history as reference category to estimate control measures' effectiveness on other parameters
    max_Inftot = max(Inftot),
    max_Inftot_time = time[which.max(Inftot)],
    # Find times at which I(t) is at a percentage of its max
    p05_Inftot_time = max(time[which(Inftot <= max_Inftot*0.05 & time < max_Inftot_time)]),
    p10_Inftot_time = max(time[which(Inftot <= max_Inftot*0.10 & time < max_Inftot_time)]),
    p25_Inftot_time = max(time[which(Inftot <= max_Inftot*0.25 & time < max_Inftot_time)]),
    p50_Inftot_time = max(time[which(Inftot <= max_Inftot*0.50 & time < max_Inftot_time)]),
    # Find times at which x number of IDX(t) are seen
    IDX500_time = max(time[which((Inftot - InfNoDX) <= 500 & time < max_Inftot_time)]),
    IDX100_time = max(time[which((Inftot - InfNoDX) <= 100 & time < max_Inftot_time)]),
    CumInfTot = sum(Inftot)) %>%
  # arrange(-ref_category_nathist) %>%
  slice_head() %>%
  ungroup() %>%
  group_by(n_hhsize, r_beta, r_tau, r_omega, n_exp_states, n_inf_states) %>%
  mutate(
    pid_nathist = max(pid*ref_category_nathist)
  ) %>%
  ungroup()

## Calculate  differential outcomes compared to reference category (Group by parameter set) ----
# Create a data.frame with only reference category (ie, hhsize = 1)
df_out_inf_all_control_summ_nathist <- df_out_inf_all_control_summ %>%
  filter(ref_category_nathist == 1) %>%
  rename(max_Inftot_nathist = max_Inftot,
         max_Inftot_time_nathist = max_Inftot_time,
         CumInfTot_nathist = CumInfTot) %>%
  dplyr::select(pid, pid_nathist, n_hhsize, n_exp_states, n_inf_states,
         r_beta, r_tau, r_omega, 
         max_Inftot_nathist, max_Inftot_time_nathist, CumInfTot_nathist)

# Compute differential outcomes compared to reference category (Group by parameter set)
df_out_inf_all_control_nathist_summ <- df_out_inf_all_control_summ %>%
  left_join(df_out_inf_all_control_summ_nathist, by = "pid_nathist") %>%
  mutate(max_Inftot_diff = max_Inftot - max_Inftot_nathist,
         max_Inftot_time_diff = max_Inftot_time - max_Inftot_time_nathist,
         CumInfTot_diff = CumInfTot - CumInfTot_nathist,
         # Percentage changes
         max_Inftot_diff_perc = (max_Inftot_diff/max_Inftot_nathist)*100,
         max_Inftot_time_diff_perc = (max_Inftot_time_diff/max_Inftot_time_nathist)*100,
         CumInfTot_diff_perc = (CumInfTot_diff/CumInfTot_nathist)*100) %>%
  rename(pid = pid.x)



## Calculate bias for all models compared bias to reference category (Group by model characteristics) ----
# Create a data.frame with only reference category (ie, hhsize = 1, e =1, i = 1)
df_out_inf_all_control_summ_ref_bias <- df_out_inf_all_control_nathist_nhh_summ %>%
  left_join(df_out_inf_all_control_nathist_summ, by = "pid") %>%
  mutate(# Absolute bias
    max_Inftot_diff_bias_abs = max_Inftot_diff_nhh - max_Inftot_diff,
    max_Inftot_time_diff_bias_abs = max_Inftot_time_diff_nhh - max_Inftot_time_diff,
    CumInfTot_diff_bias_abs = CumInfTot_diff_nhh - CumInfTot_diff,
    # Relative bias
    max_Inftot_diff_bias_rel = 100 * max_Inftot_diff_bias_abs/max_Inftot_diff,
    max_Inftot_time_diff_bias_rel = 100 * max_Inftot_time_diff_bias_abs/max_Inftot_time_diff,
    CumInfTot_diff_bias_rel = 100 * CumInfTot_diff_bias_abs/CumInfTot_diff) %>%
  rename(pid_nathist = pid_nathist.y,
         vax_prop = vax_prop.x, 
         level_npi = level_npi.x, 
         eff_vax = eff_vax.x,
         n_hhsize = n_hhsize.y.y,
         r_beta = r_beta.y.y, 
         r_tau = r_tau.y.y, 
         r_omega = r_omega.y.y,
         n_exp_states = n_exp_states.y.y,
         n_inf_states = n_inf_states.y.y,
         Esize = Esize.y,
         Isize = Isize.y,
         PropVax = PropVax.y,
         EffVax = EffVax.y,
         NPIeff = NPIeff.y,
         NPIeff_labels = NPIeff_labels.y,
         NPIeff_simple = NPIeff_simple.y,
         `Multicompartment structure` = `Multicompartment structure.y`,
         `Household size` = `Household size.y`,
         `Household size labels` = `Household size labels.y`,
         `Vaccine effectiveness` = `Vaccine effectiveness.y`) %>%
  dplyr::select(pid, pid_nathist,
         n_hhsize,
         r_beta, r_tau, r_omega, 
         vax_prop, level_npi, eff_vax,
         n_exp_states, n_inf_states,
         Esize, Isize,
         PropVax, EffVax, NPIeff,
         NPIeff_labels, NPIeff_simple,
         `Multicompartment structure`,
         `Household size`, `Household size labels`,
         `Vaccine effectiveness`,
         # CumInfTot_diff_nhh, CumInfTot_diff,
         max_Inftot_diff_bias_abs, max_Inftot_time_diff_bias_abs, CumInfTot_diff_bias_abs, 
         max_Inftot_diff_bias_rel, max_Inftot_time_diff_bias_rel, CumInfTot_diff_bias_rel)

### Control measures' effects ----
# Format data for plotting
df_control_measures_E1_I1_E3_I3_hh <- df_out_inf_all_control_summ %>% 
  filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.000 & time <= 70 & 
           ((`Multicompartment structure` %in% c("E=3, I=3") & n_hhsize %in% c(3))) & 
           eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi %in% c(0.4, 1.0)) %>%
  mutate(`Natural history structure`  = paste0(`Multicompartment structure`, " & HH = ", n_hhsize),
         NPIeff_labels = c("No NPI", "NPI60")) %>%
  dplyr::select(`Multicompartment structure`, n_hhsize, level_npi, max_Inftot, pid_nathist,
         `Natural history structure`, NPIeff, NPIeff_labels)

df_control_measures_E1_I1_E3_I3_nhh <- df_out_inf_all_control_nathist_nhh_summ %>%
  filter(pid_nathist == df_control_measures_E1_I1_E3_I3_hh$pid_nathist[1]  & 
           eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi %in% c(0.4, 1.0)) %>%
  rename(n_hhsize = n_hhsize.x) %>%
  mutate(`Natural history structure`  = paste0(`Multicompartment structure`, " & HH = ", n_hhsize),
         NPIeff_labels = c("No NPI", "NPI60")) %>%
  dplyr::select(`Multicompartment structure`, n_hhsize, level_npi, max_Inftot_nhh, pid_nathist,
         `Natural history structure`, NPIeff, NPIeff_labels) %>%
  rename(max_Inftot = max_Inftot_nhh)

df_control_measures_E1_I1_E3_I3 <- bind_rows(df_control_measures_E1_I1_E3_I3_hh,
                                             df_control_measures_E1_I1_E3_I3_nhh)

df_control_measures_E1_I1_E3_I3 

#### Figure 3: Differential effect on control measures' effects BIAS ----
png(file = "figs/Paper/Fig3_exemplary_bias.png", width = 960, height = 960*0.8)
delta_hh  <- (df_control_measures_E1_I1_E3_I3$max_Inftot[2] - df_control_measures_E1_I1_E3_I3$max_Inftot[1])
delta_nhh <- (df_control_measures_E1_I1_E3_I3$max_Inftot[4] - df_control_measures_E1_I1_E3_I3$max_Inftot[3])
rbias_fig3 <- (delta_nhh - delta_hh)/delta_hh
gg_control_measures_E1_I1_E3_I3 <- ggplot(df_control_measures_E1_I1_E3_I3,
                                          aes(x = `Natural history structure`, y = max_Inftot,
                                              group = NPIeff, fill = NPIeff_labels)) +
  geom_col(width = 0.5, position = position_dodge(0.5)) +
  # scale_fill_grey("Intervention") +
  scale_fill_manual("Intervention", values = c("No NPI" = "darkgreen", "NPI60" = "darkred")) + 
  scale_y_continuous("Peak size (thousands)", breaks = seq(0, 700000, by = 100000),
                     labels = function(x) round(x/1000, digits = 0), 
                     limits = c(0, 700000)) +
  scale_x_discrete() + 
  geom_bracket(y.position = 600*1000, label = "O[NH]", type = "expression", 
               xmin = 0.7, 
               xmax = 1.3, inherit.aes = FALSE, label.size = 8) +
  geom_bracket(y.position = 580*1000, label = "O[HH]", type = "expression", 
               xmin = 1.7, 
               xmax = 2.3, inherit.aes = FALSE, label.size = 8) +
  # xlab("Household size") +
  theme_bw(base_size = 20) +
  # coord_flip() +
  theme(strip.background = element_rect(colour = "white", fill = "white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        axis.title.x = element_blank(),
        legend.position = c(0.5, 0.3),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.key = element_blank())
gg_control_measures_E1_I1_E3_I3
## First pair of brackets
grid.brackets(#x1 = 180, x2 = 180, 
              #y1 = 60, y2 = 95,
              x1 = 430, x2 = 430, 
              y1 = 150, y2 = 285, 
              lwd = 1)
grid.text(x = unit(470, "native"), y = unit(223, "native"),
          label = expression(paste(Delta, O[NH])), hjust = 0, vjust=0, 
          gp=gpar(fontsize=20))
## Second pair of brackets
grid.brackets(x1 = 830, x2 = 830, 
              y1 = 170, y2 = 275, 
              lwd = 1)
grid.text(x = unit(870, "native"), y = unit(230, "native"),
          label = expression(paste(Delta, O[HH])), hjust = 0, vjust=0, 
          gp=gpar(fontsize=20))

grid.text(x = unit(100, "native"), y = unit(75, "native"),
          label =  expression(paste("rBias = ", 
                                    frac(paste(Delta, O[NH]) - paste(Delta, O[HH]), 
                                         paste(Delta, O[HH])
                                    )%*%100, " = 31%")
          ), 
          hjust = 0, vjust=0, 
          gp=gpar(fontsize=20))
dev.off()

#### Maximum infections ----
alphas <- c("# of I compartments = 1" = 0.4, 
            "# of I compartments = 2" = 0.7, 
            "# of I compartments = 3" = 1.0)
gg_peak_size_bias_rel <- ggplot(df_out_inf_all_control_summ_ref_bias %>% 
         filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.000 & 
                  n_hhsize %in% c(3, 5) &
                  eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi != 1), 
       aes(x =`Household size`, y = max_Inftot_diff_bias_rel, 
           group = Isize, alpha = Isize, fill = Esize)) +
  geom_bar(stat = "identity", position = position_dodge(0.8)) +
  facet_grid(NPIeff ~ Esize) +
  scale_y_continuous(breaks = seq(-20, 50, by = 5)) +
  scale_fill_discrete(l = 50) +
  scale_alpha_manual(values = alphas) +
  ylab("Relative bias on peak size (%)") +
  theme_bw(base_size = 20) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        legend.position = c(.84, .88),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.spacing.y = unit(-0.25, "cm"))
gg_peak_size_bias_rel
ggsave(plot = gg_peak_size_bias_rel, 
       filename = "figs/Paper/Fig4_control_measures_peak_size_bias_rel.pdf", 
       width = 12, height = 8)
ggsave(plot = gg_peak_size_bias_rel, 
       filename = "figs/Paper/Fig4_control_measures_peak_size_bias_rel.png", 
       width = 12, height = 8)
ggsave(plot = gg_peak_size_bias_rel, 
       filename = "figs/Paper/Fig4_control_measures_peak_size_bias_rel.jpeg", 
       width = 12, height = 8) 

my_png <- readPNG(source = "figs/Paper/Fig3_exemplary_bias.png", native = TRUE)
gg_exemplary_bias <- ggplot() + 
  background_image(raster.img = my_png)
gg_exemplary_bias

fig3 <- ggarrange(gg_exemplary_bias, gg_peak_size_bias_rel,
                  labels = c("A)", "B)"),
                  ncol = 2, nrow = 1,
                  hjust = 0, vjust = 1.5,
                  font.label = list(color = "black", size = 16))

ggsave(plot = fig3, 
       filename = "figs/Paper/Fig3_NoHH_HH.pdf", 
       width = 24, height = 10)
ggsave(plot = fig3, 
       filename = "figs/Paper/Fig3_NoHH_HH.png", 
       width = 24, height = 10)

# Appendix figures ----
## Appendix Figure 1 ----
v_fig_pid_npi20_nohh_vs_hh <- df_out_inf_all_summ %>% 
  filter(n_hhsize == 3, 
         eff_vax == 1 & 
           vax_prop == 0,
         level_npi == 0.8,
         r_beta == 0.25 & r_tau == 0.50 & 
           r_omega == 0.000 & time <= 60 &
           `Multicompartment structure` %in% c("E=1, I=1", 
                                               "E=3, I=3", 
                                               "E=3, I=1", 
                                               "E=1, I=3") 
) %>%
  dplyr::select(pid)

df_fig_npi20_nohh_vs_hh <-  bind_rows(df_out_inf_all_nhh %>% 
                                  mutate(Structure = "No HH") %>%
                                  mutate(`Household size` = ordered(1, 
                                                                    levels = c(1, 3, 5)),
                                         `NPIeff_labels` = ordered("NPI20", 
                                                                   levels = c("No NPI", "NPI20", "NPI60"))) %>%
                                  rename(#Exptot = Exptot_nhh ,
                                    InfNoDX = InfNoDX_nhh,
                                    Inftot = Inftot_nhh), 
                                df_out_inf_all %>% 
                                  mutate(Structure = "HH")) %>% 
  filter(pid %in% as.matrix(v_fig_pid_npi20_nohh_vs_hh)) %>%
  bind_rows(df_fig_nohh_vs_hh)


df_fig_npi60_20_0_nohh_vs_hh <-  bind_rows(df_out_inf_all_nhh %>% 
                                        mutate(Structure = "No HH") %>%
                                        mutate(`Household size` = ordered(1, 
                                                                          levels = c(1, 3, 5)),
                                               `NPIeff_labels` = ordered("NPI60", 
                                                                         levels = c("No NPI", "NPI20", "NPI60"))) %>%
                                        rename(#Exptot = Exptot_nhh ,
                                          InfNoDX = InfNoDX_nhh,
                                          Inftot = Inftot_nhh), 
                                      df_out_inf_all %>% 
                                        mutate(Structure = "HH")) %>% 
  filter(pid %in% as.matrix(v_fig_pid_npi60_nohh_vs_hh)) %>%
  bind_rows(df_fig_npi20_nohh_vs_hh)

gg_epidemic_curve_NPI60_20_0_E1_I1_E3_I3 <- ggplot(df_fig_npi60_20_0_nohh_vs_hh, 
                                              aes(x = time, y = Inftot/10e6, 
                                                  color = Structure, 
                                                  alpha = NPIeff_simple)) + # 
  geom_line(size = 1.1) +
  # geom_segment(aes(x = max_Inftot_time_E1_I1_hh1$max_Inftot_time, 
  #                  xend = max_Inftot_time_E1_I1_hh1$max_Inftot_time, 
  #                  y = 0, 
  #                  yend = max_Inftot_time_E1_I1_hh1$max_Inftot),
  #              linetype = "dashed", color = "blue") +
  # geom_segment(aes(x = max_Inftot_time_E1_I1_hh3$max_Inftot_time, 
  #                  xend = max_Inftot_time_E1_I1_hh3$max_Inftot_time, 
  #                  y = 0, 
  #                  yend = max_Inftot_time_E1_I1_hh3$max_Inftot),
  #              linetype = "dashed", color = "red") +
  # geom_point(aes(x = max_Inftot_time_E1_I1_hh1$max_Inftot_time, 
#                y = max_Inftot_time_E1_I1_hh1$max_Inftot),
#            color = "blue", size = 2.5) +
# geom_point(aes(x = max_Inftot_time_E1_I1_hh3$max_Inftot_time, 
#                y = max_Inftot_time_E1_I1_hh3$max_Inftot),
#            color = "red", size = 2.5) +
# geom_segment(aes(x = x0, xend = max_Inftot_time,
#                  y = y0, yend = max_Inftot,
#                  color = `Household size`),
#              linejoin = "mitre",
#              linetype = "dashed", arrow = arrow(length = unit(0.15, "inches"),
#                                                 type = "closed"), 
#              show.legend = FALSE) +
# geom_vline(xintercept = as.numeric(max_Inftot_time_E1_I1_hh1), linetype = "dashed", color = "blue") +
# geom_vline(xintercept = as.numeric(max_Inftot_time_E1_I1_hh3), linetype = "dashed", color = "red") +
facet_grid(Esize ~ Isize) +
  scale_x_continuous(limits = c(10, 50)) +
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1.0), 
                     limits = c(0, 0.06)) +
  # scale_color_grey(start = 0.2, end = 0.6) +
  scale_alpha_manual("NPI effectiveness", 
                     values = c("60%" = 1.0,
                                "20%" = 0.7, 
                                "0%"  = 0.3)) + 
  scale_color_manual(values = c("No HH" = "blue", "HH" = "red")) +
  xlab("Time") +
  ylab("Infected population (% of total population)") +
  guides(color = guide_legend(nrow = 2), 
         linetype = guide_legend(nrow = 1)) +
  theme_bw(base_size = 20) +
  theme(strip.background = element_rect(colour = "white", fill = "white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        # legend.position = c(""),
        legend.position = c(0.85, 0.80),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key = element_blank())
gg_epidemic_curve_NPI60_20_0_E1_I1_E3_I3
ggsave(plot = gg_epidemic_curve_NPI60_20_0_E1_I1_E3_I3, 
       filename = "figs/Paper/FigA1_NPI60_20_0_curves_calib.pdf", 
       width = 12, height = 8)
ggsave(plot = gg_epidemic_curve_NPI60_20_0_E1_I1_E3_I3, 
       filename = "figs/Paper/FigA1_NPI60_20_0_curves_calib.png", 
       width = 12, height = 8)

#### Appendix Figure 2 ----
v_fig_pid_vax_E1_I1_nohh_vs_hh <- df_out_inf_all_summ %>%
  filter(n_hhsize %in% c(3) & 
           eff_vax %in% c(0.5, 0.9) & 
           vax_prop %in% c(0, 0.3, 0.9) & 
           level_npi %in% c(1.0),
         r_beta == 0.25 & r_tau == 0.50 & 
           r_omega == 0.000 & time <= 60 &
           `Multicompartment structure` %in% c("E=1, I=1", 
                                               # "E=3, I=3", 
                                               # "E=3, I=1", 
                                               # "E=1, I=3",
                                               NULL) 
  ) %>%
  dplyr::select(pid)

df_epi_curve_vax_E1_I1 <- bind_rows(df_out_inf_all_nhh %>% 
                                      mutate(Structure = ordered("No HH", 
                                                                 levels = c("No HH", "HH"))) %>%
                                      mutate(`Household size` = ordered(1, 
                                                                        levels = c(1, 3, 5))) %>%
                                      rename(#Exptot = Exptot_nhh ,
                                        InfNoDX = InfNoDX_nhh,
                                        Inftot = Inftot_nhh), 
                                    df_out_inf_all %>% 
                                      mutate(Structure = ordered("HH", 
                                                                 levels = c("No HH", "HH"))) %>%
                                      filter(`Multicompartment structure` %in% c("E=1, I=1", 
                                                                                 # "E=3, I=3", 
                                                                                 # "E=3, I=1", 
                                                                                 # "E=1, I=3",
                                                                                 NULL))) %>% 
  filter(pid %in% as.matrix(v_fig_pid_vax_E1_I1_nohh_vs_hh))

gg_epidemic_curve_vax_E1_I1 <- ggplot(df_epi_curve_vax_E1_I1, 
                                      aes(x = time, y = Inftot/10e6*100, 
                                          color = as.factor(vax_prop*100))
                                      # linetype = as.factor(vax_prop))
) + # 
  geom_line(size = 1.1) +
  # geom_vline(xintercept = as.numeric(max_Inftot_time_E1_I1_hh1), linetype = "dashed", color = "blue") +
  # geom_vline(xintercept = as.numeric(max_Inftot_time_E1_I1_hh3), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "black") +
  facet_grid(EffVax ~ Structure) +
  # scale_color_grey(start = 0.2, end = 0.6) +
  # scale_color_manual(values = c("1" = "blue", "3" = "red")) +
  # scale_y_continuous(trans = "log") +
  xlim(0, 70) +
  xlab("Time") +
  ylab("Infected population (% of total population)") +
  guides(color = guide_legend(title = "Vaccine coverage (%)",
                              nrow = 3), 
         linetype = guide_legend(title = "Vaccine coverage (%)",
                                 nrow = 1)) +
  theme_bw(base_size = 20) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        # legend.position = c(""),
        legend.position = c(0.85, 0.90),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key = element_blank())
gg_epidemic_curve_vax_E1_I1
ggsave(plot = gg_epidemic_curve_vax_E1_I1, 
       filename = "figs/Paper/FigA3_epi_curves_vax_E1_I1_calib.pdf", 
       width = 12, height = 8)

#### Appendix Figure 3 ----
v_fig_pid_vax_E3_I3_nohh_vs_hh <- df_out_inf_all_summ %>%
  filter(n_hhsize %in% c(1, 3) & 
           eff_vax %in% c(0.5, 0.9) & 
           vax_prop %in% c(0, 0.3, 0.9) & 
           level_npi %in% c(1.0),
         r_beta == 0.25 & r_tau == 0.50 & 
           r_omega == 0.000 & time <= 60 &
           `Multicompartment structure` %in% c(#"E=1, I=1", 
                                               "E=3, I=3",
                                               # "E=3, I=1", 
                                               # "E=1, I=3",
                                               NULL) 
  ) %>%
  dplyr::select(pid)

df_epi_curve_vax_E3_I3 <- bind_rows(df_out_inf_all_nhh %>% 
                                      mutate(Structure = ordered("No HH", 
                                                                 levels = c("No HH", "HH"))) %>%
                                      mutate(`Household size` = ordered(1, 
                                                                        levels = c(1, 3, 5))) %>%
                                      rename(#Exptot = Exptot_nhh ,
                                        InfNoDX = InfNoDX_nhh,
                                        Inftot = Inftot_nhh), 
                                    df_out_inf_all %>% 
                                      mutate(Structure = ordered("HH", 
                                                                 levels = c("No HH", "HH"))) %>%
                                      filter(n_hhsize == 3, 
                                             `Multicompartment structure` %in% c(#"E=1, I=1", 
                                                                                 "E=3, I=3", 
                                                                                 # "E=3, I=1", 
                                                                                 # "E=1, I=3",
                                                                                 NULL))) %>% 
  filter(pid %in% as.matrix(v_fig_pid_vax_E3_I3_nohh_vs_hh))

gg_epidemic_curve_vax_E3_I3 <- ggplot(df_epi_curve_vax_E3_I3, 
                                      aes(x = time, y = Inftot/10e6*100, 
                                          color = as.factor(vax_prop*100))
                                      # linetype = as.factor(vax_prop))
) + # 
  geom_line(size = 1.1) +
  # geom_vline(xintercept = as.numeric(max_Inftot_time_E1_I1_hh1), linetype = "dashed", color = "blue") +
  # geom_vline(xintercept = as.numeric(max_Inftot_time_E1_I1_hh3), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "black") +
  facet_grid(EffVax ~ Structure) +
  # scale_color_grey(start = 0.2, end = 0.6) +
  # scale_color_manual(values = c("1" = "blue", "3" = "red")) +
  # scale_y_continuous(trans = "log") +
  xlim(0, 70) +
  xlab("Time") +
  ylab("Infected population (% of total population)") +
  guides(color = guide_legend(title = "Vaccine coverage (%)",
                              nrow = 3), 
         linetype = guide_legend(title = "Vaccine coverage (%)",
                                 nrow = 1)) +
  theme_bw(base_size = 20) +
  theme(strip.background = element_rect(colour="white", fill = "white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        # legend.position = c(""),
        legend.position = c(0.85, 0.90),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key = element_blank())
gg_epidemic_curve_vax_E3_I3
ggsave(plot = gg_epidemic_curve_vax_E3_I3, 
       filename = "figs/Paper/FigA4_epi_curves_vax_E3_I3_calib.pdf", 
       width = 12, height = 8)

#### Time of epidemic peak infections ----
ggplot(df_out_inf_all_control_summ_ref_bias %>% 
         filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.000 & 
                  n_hhsize %in% c(3, 5) &
                  eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi != 1), 
       aes(x =`Household size`, y = max_Inftot_time_diff_bias_rel, 
           group = Isize, alpha = Isize, fill = Esize)) +
  geom_bar(stat = "identity", position = position_dodge(0.8)) +
  facet_grid(NPIeff ~ Esize) +
  scale_fill_discrete(l = 50) +
  scale_alpha_manual(values = alphas) +
  ylab("Relative bias on peak time (%)") +
  theme_bw(base_size = 20) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        legend.position = c(.84, .86),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank())

#### Epidemic Size ----
ggplot(df_out_inf_all_control_summ_ref_bias %>% 
         filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.000 & 
                  n_hhsize %in% c(3, 5) &
                  eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi != 1), 
       aes(x =`Household size`, y = CumInfTot_diff_bias_rel, 
           group = Isize, alpha = Isize, fill = Esize)) +
  geom_bar(stat = "identity", position = position_dodge(0.8)) +
  facet_grid(NPIeff ~ Esize) +
  scale_fill_discrete(l = 50) +
  scale_alpha_manual(values = alphas) +
  ylab("Relative bias on epidemic size (%)") +
  theme_bw(base_size = 20) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        legend.position = c(.84, .86),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank())

#### Appendix Figure 4 ----
# Format data for plotting
df_control_measures_bias_vax <- df_out_inf_all_control_summ_ref_bias %>% 
  filter(r_beta == 0.25 & r_tau == 0.50 & 
           r_omega == 0.000 & 
           `Multicompartment structure` %in% c("E=1, I=1",
                                               "E=3, I=3",
                                               "E=3, I=1",
                                               "E=1, I=3"
           ) & 
           n_hhsize %in% c(1, 3, 5) & 
           eff_vax %in% c(0.5, 0.9) & 
           vax_prop %in% c(0.3, 0.9) & 
           level_npi %in% c(1.0))

# df_control_measures_bias_vax$CumInfTot_diff_nhh - df_control_measures_bias_vax$CumInfTot_diff
# df_control_measures_bias_vax$CumInfTot_diff_bias_abs
# df_control_measures_bias_vax$CumInfTot_diff_bias_abs/df_control_measures_bias_vax$CumInfTot_diff
# hist(df_control_measures_bias_vax$CumInfTot_diff_bias_rel)

gg_control_measures_bias_vax_CumInfTot_diff_bias_rel <- ggplot(df_control_measures_bias_vax, 
                                                               aes(x = `Household size`, 
                                                                   fill = PropVax,
                                                                   y = CumInfTot_diff_bias_rel)) +
  geom_col(width = 0.5, position = position_dodge(0.5)) +
  facet_grid(`EffVax` ~ `Multicompartment structure`) +
  scale_y_continuous("Relative bias (%)") +
  theme_bw(base_size = 20) +
  theme(legend.position = "bottom", 
        legend.title = element_blank())
gg_control_measures_bias_vax_CumInfTot_diff_bias_rel
ggsave(plot = gg_control_measures_bias_vax_CumInfTot_diff_bias_rel, 
       filename = "figs/Paper/FigA5_control_measures_bias_vax_CumInfTot_diff_bias_rel_calib.pdf", 
       width = 12, height = 9)
ggsave(plot = gg_control_measures_bias_vax_CumInfTot_diff_bias_rel, 
       filename = "figs/Paper/FigA5_control_measures_bias_vax_CumInfTot_diff_bias_rel_calib.png", 
       width = 12, height = 9)

#### Appendix Figure 5 ----
gg_control_measures_bias_vax_MaxInfTot_diff_bias_rel <- ggplot(df_control_measures_bias_vax, 
                                                               aes(x = `Household size`, 
                                                                   fill = PropVax,
                                                                   y = max_Inftot_diff_bias_rel)) +
  geom_col(width = 0.5, position = position_dodge(0.5)) +
  facet_grid(`EffVax` ~ `Multicompartment structure`) +
  scale_y_continuous(n.breaks = 6, "Relative bias (%)") +
  theme_bw(base_size = 20) +
  theme(legend.position = "bottom", 
        legend.title = element_blank())
gg_control_measures_bias_vax_MaxInfTot_diff_bias_rel
ggsave(plot = gg_control_measures_bias_vax_MaxInfTot_diff_bias_rel, 
       filename = "figs/Paper/FigA6_control_measures_bias_vax_MaxInfTot_diff_bias_rel_calib.pdf", 
       width = 12, height = 9)
ggsave(plot = gg_control_measures_bias_vax_MaxInfTot_diff_bias_rel, 
       filename = "figs/Paper/FigA6_control_measures_bias_vax_MaxInfTot_diff_bias_rel_calib.png", 
       width = 12, height = 9)

### Table 3: Meta regression on category-specific outcomes on Control Measures NPI 60% reduction ----
df_out_inf_all_control_summ_ref_bias_metareg <- df_out_inf_all_control_summ_ref_bias %>% 
  filter(r_omega == 0.000 & 
           !((`Multicompartment structure` %in% c("E=1, I=1") & n_hhsize %in% c(1))) & 
           eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi %in% c(0.4)) 
df_out_inf_all_control_summ_ref_bias_metareg$HH_5 <- df_out_inf_all_control_summ_ref_bias_metareg$n_hhsize == 5

## Without interactions
fit_hh_peak_time_bias_rel <- lm(max_Inftot_time_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                  HH_5 +
                                  r_tau + r_beta, 
                                data = df_out_inf_all_control_summ_ref_bias_metareg %>% 
                                  filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel)

fit_hh_peak_size_bias_rel <- lm(max_Inftot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                  HH_5 +
                                  r_tau + r_beta,
                                data = df_out_inf_all_control_summ_ref_bias_metareg)
summary(fit_hh_peak_size_bias_rel)

fit_hh_epidemic_size_rel <- lm(CumInfTot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                 HH_5 +
                                 r_tau + r_beta, 
                               data = df_out_inf_all_control_summ_ref_bias_metareg)
summary(fit_hh_epidemic_size_rel)

# Linear combination coefficient tests
l_lincomm_rel <- vector(mode = "list", length = 3)
l_lincomm_rel[[1]] <- summary(multcomp::glht(fit_hh_peak_time_bias_rel, linfct = c("(Intercept) + HH_5TRUE = 0")))
l_lincomm_rel[[2]] <- summary(multcomp::glht(fit_hh_peak_size_bias_rel, linfct = c("(Intercept) + HH_5TRUE = 0")))
l_lincomm_rel[[3]] <- summary(multcomp::glht(fit_hh_epidemic_size_rel, linfct = c("(Intercept) + HH_5TRUE = 0")))

df_lincomm_rel <- data.frame(`Epidemic outcome` = c("Peak time", "Peak size", "Epidemic size"),
                                 coef = as.numeric(c(l_lincomm_rel[[1]]$test$coefficients, 
                                                     l_lincomm_rel[[2]]$test$coefficients,
                                                     l_lincomm_rel[[3]]$test$coefficients)), 
                                 pval = c(l_lincomm_rel[[1]]$test$pvalues[1], 
                                          l_lincomm_rel[[2]]$test$pvalues[1],
                                          l_lincomm_rel[[3]]$test$pvalues[1]), 
                                 check.names = F)
df_lincomm_rel
xtable::xtable(df_lincomm_rel)

## With interactions
fit_hh_peak_time_bias_rel_int <- lm(max_Inftot_time_diff_bias_rel ~ #n_exp_states*n_inf_states + 
                                      n_exp_states*HH_5 + 
                                      n_inf_states*HH_5 +
                                      r_tau + r_beta, 
                                    data = df_out_inf_all_control_summ_ref_bias_metareg %>% 
                                      filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel_int)

fit_hh_peak_size_bias_rel_int <- lm(max_Inftot_diff_bias_rel ~ #n_exp_states*n_inf_states + 
                                      n_exp_states*HH_5 + 
                                      n_inf_states*HH_5 +
                                      r_tau + r_beta,
                                    data = df_out_inf_all_control_summ_ref_bias_metareg)
summary(fit_hh_peak_size_bias_rel_int)

fit_hh_epidemic_size_rel_int <- lm(CumInfTot_diff_bias_rel ~ #n_exp_states*n_inf_states + 
                                     n_exp_states*HH_5 + 
                                     n_inf_states*HH_5 +
                                     r_tau + r_beta, 
                                   data = df_out_inf_all_control_summ_ref_bias_metareg)
summary(fit_hh_epidemic_size_rel_int)

# Linear combination coefficient tests
l_lincomm_rel_int <- vector(mode = "list", length = 3)
l_lincomm_rel_int[[1]] <- summary(multcomp::glht(fit_hh_peak_time_bias_rel_int, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))
l_lincomm_rel_int[[2]] <- summary(multcomp::glht(fit_hh_peak_size_bias_rel_int, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))
l_lincomm_rel_int[[3]] <- summary(multcomp::glht(fit_hh_epidemic_size_rel_int, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))

df_lincomm_rel_int <- data.frame(`Epidemic outcome` = c("Peak time", "Peak size", "Epidemic size"),
                                 coef = as.numeric(c(l_lincomm_rel_int[[1]]$test$coefficients, 
                                                     l_lincomm_rel_int[[2]]$test$coefficients,
                                                     l_lincomm_rel_int[[3]]$test$coefficients)), 
                                 pval = c(l_lincomm_rel_int[[1]]$test$pvalues[1], 
                                          l_lincomm_rel_int[[2]]$test$pvalues[1],
                                          l_lincomm_rel_int[[3]]$test$pvalues[1]), 
                                 check.names = F)
df_lincomm_rel_int
xtable::xtable(df_lincomm_rel_int)

df_lincomm_rel$Interaction <- "No"
df_lincomm_rel_int$Interaction <- "Yes"

df_lincomm_rel_comb <- bind_rows(df_lincomm_rel, df_lincomm_rel_int) %>%
  group_by(`Epidemic outcome`) %>%
  arrange(desc(`Epidemic outcome`), Interaction)
df_lincomm_rel_nhh_comb

v_pval_rel <- df_lincomm_rel_comb$pval
v_star_rel <- cut(v_pval_rel,
                      breaks = c(-Inf, 0.01, 0.05, 0.1, Inf), 
                      labels = c("$^{***}$", "$^{**}$", "$^{*}$", ""))

v_coef_rel <- as.character(round(df_lincomm_rel_comb$coef, 3))
v_coef_sg_rel <- paste0(v_coef_rel, v_star_rel)

v_stargazer_int_effects_order_main <- c(8, 3, 1, 2, 5, 4, 6, 7)
v_stargazer_int_effects_label_main <- c("Household size 3 (HH3)", "Household size 5 (HH5)",
                                        "E", "I", "$\\tau$", "$\\beta$",
                                        "E*HH5", "I*HH5")
v_stargazer_notes_int_effects <- c( "HH3 describes how much the intervention effect differs if household structure is excluded",
                       "from an otherwise similar epidemic where the true household size is 3.",
                       "HH5 describes the incremental difference in the intervention effect if household structure is excluded",
                       "from an otherwise similar epidemic where the true household size is 5 instead of 3.",
                       "For the models, the total effect of exclusion of household structure involves interactions with other terms.",
                       "For example, E*HH5 describes how the incremental difference depends upon the number of exposed",
                       "compartments when the true household size is 5 instead of 3.",
                       "We estimated the magnitude of the linear combination of the relevant coefficients and tested their",
                       "significance.")

stargazer(fit_hh_peak_time_bias_rel, fit_hh_peak_time_bias_rel_int,
          fit_hh_peak_size_bias_rel, fit_hh_peak_size_bias_rel_int,
          fit_hh_epidemic_size_rel, fit_hh_epidemic_size_rel_int,
          type = "text", 
          title = "Metaregression estimates on relative bias of treatment effects, NPI = 60%",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          label = "tab:metaregbias",
          order = v_stargazer_int_effects_order_main,
          covariate.labels = v_stargazer_int_effects_label_main,
          notes = v_stargazer_notes_int_effects,
          add.lines = list(c("HH3 + HH5", 
                             v_coef_sg_rel)),
          notes.append = TRUE,
          notes.align = "l",
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE)

stargazer(fit_hh_peak_time_bias_rel, fit_hh_peak_time_bias_rel_int,
          fit_hh_peak_size_bias_rel, fit_hh_peak_size_bias_rel_int,
          fit_hh_epidemic_size_rel, fit_hh_epidemic_size_rel_int,
          type = "latex", 
          out = "output/Table3_metaregression_bias.tex", 
          title = "Metaregression estimates on relative bias of treatment effects, NPI = 60%",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          label = "tab:metaregbias",
          order = v_stargazer_int_effects_order_main,
          covariate.labels = v_stargazer_int_effects_label_main,
          notes = v_stargazer_notes_int_effects,
          add.lines = list(c("HH3 + HH5", 
                             v_coef_sg_rel)),
          notes.append = TRUE,
          notes.align = "l",
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE)


### Table A1: Meta regression on category-specific outcomes on Control Measures NPI 20% reduction ----
df_out_inf_all_control_summ_ref_bias$HH_5 <- df_out_inf_all_control_summ_ref_bias$n_hhsize == 5

df_out_inf_all_control_summ_ref_bias_metareg_NPI20 <- df_out_inf_all_control_summ_ref_bias %>% 
  filter(r_omega == 0.000 & 
           !((`Multicompartment structure` %in% c("E=1, I=1") & n_hhsize %in% c(1))) & 
           eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi %in% c(0.8)) 

## Without interactions
fit_hh_peak_time_bias_rel_NPI20 <- lm(max_Inftot_time_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                    HH_5 +
                                  r_tau + r_beta, 
                                data = df_out_inf_all_control_summ_ref_bias_metareg_NPI20 %>% 
                                  filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel_NPI20)

fit_hh_peak_size_bias_rel_NPI20 <- lm(max_Inftot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                  HH_5 +
                                  r_tau + r_beta,
                                data = df_out_inf_all_control_summ_ref_bias_metareg_NPI20)
summary(fit_hh_peak_size_bias_rel_NPI20)

fit_hh_epidemic_size_rel_NPI20 <- lm(CumInfTot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                       HH_5 +
                                   r_tau + r_beta, 
                               data = df_out_inf_all_control_summ_ref_bias_metareg_NPI20)
summary(fit_hh_epidemic_size_rel_NPI20)
# Linear combination coefficient tests
l_lincomm_rel_NPI20 <- vector(mode = "list", length = 3)
l_lincomm_rel_NPI20[[1]] <- summary(multcomp::glht(fit_hh_peak_time_bias_rel_NPI20, linfct = c("(Intercept) + HH_5TRUE = 0")))
l_lincomm_rel_NPI20[[2]] <- summary(multcomp::glht(fit_hh_peak_size_bias_rel_NPI20, linfct = c("(Intercept) + HH_5TRUE = 0")))
l_lincomm_rel_NPI20[[3]] <- summary(multcomp::glht(fit_hh_epidemic_size_rel_NPI20, linfct = c("(Intercept) + HH_5TRUE = 0")))

df_lincomm_rel_NPI20 <- data.frame(`Epidemic outcome` = c("Peak time", "Peak size", "Epidemic size"),
                             coef = as.numeric(c(l_lincomm_rel_NPI20[[1]]$test$coefficients, 
                                                 l_lincomm_rel_NPI20[[2]]$test$coefficients,
                                                 l_lincomm_rel_NPI20[[3]]$test$coefficients)), 
                             pval = c(l_lincomm_rel_NPI20[[1]]$test$pvalues[1], 
                                      l_lincomm_rel_NPI20[[2]]$test$pvalues[1],
                                      l_lincomm_rel_NPI20[[3]]$test$pvalues[1]), 
                             check.names = F)
df_lincomm_rel_NPI20
xtable::xtable(df_lincomm_rel_NPI20)

## With interactions
fit_hh_peak_time_bias_rel_int_NPI20 <- lm(max_Inftot_time_diff_bias_rel ~ # n_exp_states*n_inf_states + 
                                            n_exp_states*HH_5 + 
                                            n_inf_states*HH_5 +
                                      r_tau + r_beta, 
                                    data = df_out_inf_all_control_summ_ref_bias_metareg_NPI20 %>% 
                                      filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel_int_NPI20)

fit_hh_peak_size_bias_rel_int_NPI20 <- lm(max_Inftot_diff_bias_rel ~ # n_exp_states*n_inf_states + 
                                            n_exp_states*HH_5 + 
                                            n_inf_states*HH_5 +
                                      r_tau + r_beta,
                                    data = df_out_inf_all_control_summ_ref_bias_metareg_NPI20)
summary(fit_hh_peak_size_bias_rel_int_NPI20)

fit_hh_epidemic_size_rel_int_NPI20 <- lm(CumInfTot_diff_bias_rel ~ # n_exp_states*n_inf_states + 
                                           n_exp_states*HH_5 + 
                                           n_inf_states*HH_5 +
                                     r_tau + r_beta, 
                                   data = df_out_inf_all_control_summ_ref_bias_metareg_NPI20)
summary(fit_hh_epidemic_size_rel_int_NPI20)

# Linear combination coefficient tests
l_lincomm_rel_int_NPI20 <- vector(mode = "list", length = 3)
l_lincomm_rel_int_NPI20[[1]] <- summary(multcomp::glht(fit_hh_peak_time_bias_rel_int_NPI20, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))
l_lincomm_rel_int_NPI20[[2]] <- summary(multcomp::glht(fit_hh_peak_size_bias_rel_int_NPI20, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))
l_lincomm_rel_int_NPI20[[3]] <- summary(multcomp::glht(fit_hh_epidemic_size_rel_int_NPI20, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))

df_lincomm_rel_int_NPI20 <- data.frame(`Epidemic outcome` = c("Peak time", "Peak size", "Epidemic size"),
                             coef = as.numeric(c(l_lincomm_rel_int_NPI20[[1]]$test$coefficients, 
                                                 l_lincomm_rel_int_NPI20[[2]]$test$coefficients,
                                                 l_lincomm_rel_int_NPI20[[3]]$test$coefficients)), 
                             pval = c(l_lincomm_rel_int_NPI20[[1]]$test$pvalues[1], 
                                      l_lincomm_rel_int_NPI20[[2]]$test$pvalues[1],
                                      l_lincomm_rel_int_NPI20[[3]]$test$pvalues[1]), 
                             check.names = F)
df_lincomm_rel_int_NPI20
xtable::xtable(df_lincomm_rel_int_NPI20)

df_lincomm_rel_NPI20$Interaction <- "No"
df_lincomm_rel_int_NPI20$Interaction <- "Yes"

df_lincomm_rel_NPI20_comb <- bind_rows(df_lincomm_rel_NPI20, 
                                       df_lincomm_rel_int_NPI20) %>%
  group_by(`Epidemic outcome`) %>%
  arrange(desc(`Epidemic outcome`), Interaction)
df_lincomm_rel_NPI20_comb

v_pval_rel_NPI20 <- df_lincomm_rel_NPI20_comb$pval
v_star_rel_NPI20 <- cut(v_pval_rel_NPI20,
                  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf), 
                  labels = c("$^{***}$", "$^{**}$", "$^{*}$", ""))

v_coef_rel_NPI20 <- as.character(round(df_lincomm_rel_NPI20_comb$coef, 3))
v_coef_sg_rel_NPI20 <- paste0(v_coef_rel_NPI20, v_star_rel_NPI20)

stargazer(fit_hh_peak_time_bias_rel_NPI20, fit_hh_peak_time_bias_rel_int_NPI20,
          fit_hh_peak_size_bias_rel_NPI20, fit_hh_peak_size_bias_rel_int_NPI20,
          fit_hh_epidemic_size_rel_NPI20, fit_hh_epidemic_size_rel_int_NPI20,
          type = "text", 
          title = "Metaregression estimates on relative bias of treatment effects, NPI = 20%",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          label = "tab:metaregbias_app",
          order = v_stargazer_int_effects_order_main,
          covariate.labels = v_stargazer_int_effects_label_main,
          notes = v_stargazer_notes_int_effects,
          add.lines = list(c("HH3 + HH5", 
                             v_coef_sg_rel_NPI20)),
          notes.append = TRUE,
          notes.align = "l",
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE)

stargazer(fit_hh_peak_time_bias_rel_NPI20, fit_hh_peak_time_bias_rel_int_NPI20,
          fit_hh_peak_size_bias_rel_NPI20, fit_hh_peak_size_bias_rel_int_NPI20,
          fit_hh_epidemic_size_rel_NPI20, fit_hh_epidemic_size_rel_int_NPI20,
          type = "latex", 
          out = "output/TableA1_metaregression_bias_NPI20.tex", 
          title = "Metaregression estimates on relative bias of treatment effects, NPI = 20%",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          label = "tab:metaregbias_app",
          order = v_stargazer_int_effects_order_main,
          covariate.labels = v_stargazer_int_effects_label_main,
          notes = v_stargazer_notes_int_effects,
          add.lines = list(c("HH3 + HH5", 
                             v_coef_sg_rel_NPI20)),
          notes.append = TRUE,
          notes.align = "l",
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE)

### Table A2: Meta regression on category-specific outcomes on Vaccine Measures Coverage = 30% and Effectiveness = 90% NPI 0% reduction ----
df_control_measures_bias_vax_propL_effH <- df_out_inf_all_control_summ_ref_bias %>% 
  filter(r_omega == 0.000 & 
           !((`Multicompartment structure` %in% c("E=1, I=1") & n_hhsize %in% c(1))) & 
           eff_vax %in% c(0.9) & vax_prop == 0.3, level_npi %in% c(1.0)) 

## Without interactions
fit_hh_peak_time_bias_rel_vax_propL_effH <- lm(max_Inftot_time_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                                 HH_5 + 
                                                 r_tau + r_beta, 
                                               data = df_control_measures_bias_vax_propL_effH %>% 
                                                 filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel_vax_propL_effH)

fit_hh_peak_size_bias_rel_vax_propL_effH <- lm(max_Inftot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                                 HH_5 + 
                                                 r_tau + r_beta,
                                               data = df_control_measures_bias_vax_propL_effH)
summary(fit_hh_peak_size_bias_rel_vax_propL_effH)

fit_hh_epidemic_size_rel_vax_propL_effH <- lm(CumInfTot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                                HH_5 + 
                                                r_tau + r_beta, 
                                              data = df_control_measures_bias_vax_propL_effH)
summary(fit_hh_epidemic_size_rel_vax_propL_effH)
# Linear combination coefficient tests
l_lincomm_rel_vax_propL_effH <- vector(mode = "list", length = 3)
l_lincomm_rel_vax_propL_effH[[1]] <- summary(multcomp::glht(fit_hh_peak_time_bias_rel_vax_propL_effH, linfct = c("(Intercept) + HH_5TRUE = 0")))
l_lincomm_rel_vax_propL_effH[[2]] <- summary(multcomp::glht(fit_hh_peak_size_bias_rel_vax_propL_effH, linfct = c("(Intercept) + HH_5TRUE = 0")))
l_lincomm_rel_vax_propL_effH[[3]] <- summary(multcomp::glht(fit_hh_epidemic_size_rel_vax_propL_effH, linfct = c("(Intercept) + HH_5TRUE = 0")))

df_lincomm_rel_vax_propL_effH <- data.frame(`Epidemic outcome` = c("Peak time", "Peak size", "Epidemic size"),
                                   coef = as.numeric(c(l_lincomm_rel_vax_propL_effH[[1]]$test$coefficients, 
                                                       l_lincomm_rel_vax_propL_effH[[2]]$test$coefficients,
                                                       l_lincomm_rel_vax_propL_effH[[3]]$test$coefficients)), 
                                   pval = c(l_lincomm_rel_vax_propL_effH[[1]]$test$pvalues[1], 
                                            l_lincomm_rel_vax_propL_effH[[2]]$test$pvalues[1],
                                            l_lincomm_rel_vax_propL_effH[[3]]$test$pvalues[1]), 
                                   check.names = F)
df_lincomm_rel_vax_propL_effH
xtable::xtable(df_lincomm_rel_vax_propL_effH)

## With interactions
fit_hh_peak_time_bias_rel_int_vax_propL_effH <- lm(max_Inftot_time_diff_bias_rel ~ # n_exp_states*n_inf_states + 
                                                     n_exp_states*HH_5 + 
                                                     n_inf_states*HH_5 + 
                                                     r_tau + r_beta, 
                                                   data = df_control_measures_bias_vax_propL_effH %>% 
                                                     filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel_int_vax_propL_effH)

fit_hh_peak_size_bias_rel_int_vax_propL_effH <- lm(max_Inftot_diff_bias_rel ~ # n_exp_states*n_inf_states + 
                                                     n_exp_states*HH_5 + 
                                                     n_inf_states*HH_5 + 
                                                     r_tau + r_beta,
                                                   data = df_control_measures_bias_vax_propL_effH)
summary(fit_hh_peak_size_bias_rel_int_vax_propL_effH)

fit_hh_epidemic_size_rel_int_vax_propL_effH <- lm(CumInfTot_diff_bias_rel ~ # n_exp_states*n_inf_states + 
                                                    n_exp_states*HH_5 +  
                                                    n_inf_states*HH_5 +  
                                                    r_tau + r_beta, 
                                                  data = df_control_measures_bias_vax_propL_effH)
summary(fit_hh_epidemic_size_rel_int_vax_propL_effH)
# Linear combination coefficient tests
l_lincomm_rel_int_vax_propL_effH <- vector(mode = "list", length = 3)
l_lincomm_rel_int_vax_propL_effH[[1]] <- summary(multcomp::glht(fit_hh_peak_time_bias_rel_int_vax_propL_effH, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))
l_lincomm_rel_int_vax_propL_effH[[2]] <- summary(multcomp::glht(fit_hh_peak_size_bias_rel_int_vax_propL_effH, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))
l_lincomm_rel_int_vax_propL_effH[[3]] <- summary(multcomp::glht(fit_hh_epidemic_size_rel_int_vax_propL_effH, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))

df_lincomm_rel_int_vax_propL_effH <- data.frame(`Epidemic outcome` = c("Peak time", "Peak size", "Epidemic size"),
                                       coef = as.numeric(c(l_lincomm_rel_int_vax_propL_effH[[1]]$test$coefficients, 
                                                           l_lincomm_rel_int_vax_propL_effH[[2]]$test$coefficients,
                                                           l_lincomm_rel_int_vax_propL_effH[[3]]$test$coefficients)), 
                                       pval = c(l_lincomm_rel_int_vax_propL_effH[[1]]$test$pvalues[1], 
                                                l_lincomm_rel_int_vax_propL_effH[[2]]$test$pvalues[1],
                                                l_lincomm_rel_int_vax_propL_effH[[3]]$test$pvalues[1]), 
                                       check.names = F)
df_lincomm_rel_int_vax_propL_effH
xtable::xtable(df_lincomm_rel_int_vax_propL_effH)

df_lincomm_rel_vax_propL_effH$Interaction <- "No"
df_lincomm_rel_int_vax_propL_effH$Interaction <- "Yes"

df_lincomm_rel_vax_propL_effH_comb <- bind_rows(df_lincomm_rel_vax_propL_effH, 
                                       df_lincomm_rel_int_vax_propL_effH) %>%
  group_by(`Epidemic outcome`) %>%
  arrange(desc(`Epidemic outcome`), Interaction)
df_lincomm_rel_vax_propL_effH_comb

v_pval_rel_vax_propL_effH <- df_lincomm_rel_vax_propL_effH_comb$pval
v_star_rel_vax_propL_effH <- cut(v_pval_rel_vax_propL_effH,
                        breaks = c(-Inf, 0.01, 0.05, 0.1, Inf), 
                        labels = c("$^{***}$", "$^{**}$", "$^{*}$", ""))

v_coef_rel_vax_propL_effH <- as.character(round(df_lincomm_rel_vax_propL_effH_comb$coef, 3))
v_coef_sg_rel_vax_propL_effH <- paste0(v_coef_rel_vax_propL_effH, v_star_rel_vax_propL_effH)

stargazer(fit_hh_peak_time_bias_rel_vax_propL_effH, fit_hh_peak_time_bias_rel_int_vax_propL_effH,
          fit_hh_peak_size_bias_rel_vax_propL_effH, fit_hh_peak_size_bias_rel_int_vax_propL_effH,
          fit_hh_epidemic_size_rel_vax_propL_effH, fit_hh_epidemic_size_rel_int_vax_propL_effH,
          type = "text", 
          title = "Metaregression estimates vaccination coverage = 30% & effectiveness = 90%",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          label = "tab:vax_metareg_30_90",
          order = v_stargazer_int_effects_order_main,
          covariate.labels = v_stargazer_int_effects_label_main,
          add.lines = list(c("HH3 + HH5", 
                             v_coef_sg_rel_vax_propL_effH)),
          notes = v_stargazer_notes_int_effects,
          notes.append = TRUE,
          notes.align = "l",
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE)

stargazer(fit_hh_peak_time_bias_rel_vax_propL_effH, fit_hh_peak_time_bias_rel_int_vax_propL_effH,
          fit_hh_peak_size_bias_rel_vax_propL_effH, fit_hh_peak_size_bias_rel_int_vax_propL_effH,
          fit_hh_epidemic_size_rel_vax_propL_effH, fit_hh_epidemic_size_rel_int_vax_propL_effH,
          type = "latex", 
          out = "output/TableA2_metaregression_bias_vax_propL_effH.tex", 
          title = "Metaregression estimates vaccination coverage = 30% & effectiveness = 90%",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          label = "tab:vax_metareg_30_90",
          order = v_stargazer_int_effects_order_main,
          covariate.labels = v_stargazer_int_effects_label_main,
          add.lines = list(c("HH3 + HH5", 
                             v_coef_sg_rel_vax_propL_effH)),
          notes = v_stargazer_notes_int_effects,
          notes.append = TRUE,
          notes.align = "l",
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE)

### Table A3: Meta regression on category-specific outcomes on Vaccine Measures Coverage = 90% and Effectiveness = 90% NPI 0% reduction ----
df_control_measures_bias_vax_propH_effH <- df_out_inf_all_control_summ_ref_bias %>% 
  filter(r_omega == 0.000 & 
           !((`Multicompartment structure` %in% c("E=1, I=1") & n_hhsize %in% c(1))) & 
           eff_vax %in% c(0.9) & vax_prop == 0.9, level_npi %in% c(1.0)) 

## Without interactions
fit_hh_peak_time_bias_rel_vax_propH_effH <- lm(max_Inftot_time_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                                 HH_5 + 
                                                 r_tau + r_beta, 
                                               data = df_control_measures_bias_vax_propH_effH %>% 
                                                 filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel_vax_propH_effH)

fit_hh_peak_size_bias_rel_vax_propH_effH <- lm(max_Inftot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                                 HH_5 + 
                                                 r_tau + r_beta,
                                               data = df_control_measures_bias_vax_propH_effH)
summary(fit_hh_peak_size_bias_rel_vax_propH_effH)

fit_hh_epidemic_size_rel_vax_propH_effH <- lm(CumInfTot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                                HH_5 + 
                                                r_tau + r_beta, 
                                              data = df_control_measures_bias_vax_propH_effH)
summary(fit_hh_epidemic_size_rel_vax_propH_effH)
# Linear combination coefficient tests
l_lincomm_rel_vax_propH_effH <- vector(mode = "list", length = 3)
l_lincomm_rel_vax_propH_effH[[1]] <- summary(multcomp::glht(fit_hh_peak_time_bias_rel_vax_propH_effH, linfct = c("(Intercept) + HH_5TRUE = 0")))
l_lincomm_rel_vax_propH_effH[[2]] <- summary(multcomp::glht(fit_hh_peak_size_bias_rel_vax_propH_effH, linfct = c("(Intercept) + HH_5TRUE = 0")))
l_lincomm_rel_vax_propH_effH[[3]] <- summary(multcomp::glht(fit_hh_epidemic_size_rel_vax_propH_effH, linfct = c("(Intercept) + HH_5TRUE = 0")))

df_lincomm_rel_vax_propH_effH <- data.frame(`Epidemic outcome` = c("Peak time", "Peak size", "Epidemic size"),
                                            coef = as.numeric(c(l_lincomm_rel_vax_propH_effH[[1]]$test$coefficients, 
                                                                l_lincomm_rel_vax_propH_effH[[2]]$test$coefficients,
                                                                l_lincomm_rel_vax_propH_effH[[3]]$test$coefficients)), 
                                            pval = c(l_lincomm_rel_vax_propH_effH[[1]]$test$pvalues[1], 
                                                     l_lincomm_rel_vax_propH_effH[[2]]$test$pvalues[1],
                                                     l_lincomm_rel_vax_propH_effH[[3]]$test$pvalues[1]), 
                                            check.names = F)
df_lincomm_rel_vax_propH_effH
xtable::xtable(df_lincomm_rel_vax_propH_effH)

## With interactions
fit_hh_peak_time_bias_rel_int_vax_propH_effH <- lm(max_Inftot_time_diff_bias_rel ~ # n_exp_states*n_inf_states + 
                                                     n_exp_states*HH_5 + 
                                                     n_inf_states*HH_5 + 
                                                     r_tau + r_beta, 
                                                   data = df_control_measures_bias_vax_propH_effH %>% 
                                                     filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel_int_vax_propH_effH)

fit_hh_peak_size_bias_rel_int_vax_propH_effH <- lm(max_Inftot_diff_bias_rel ~ # n_exp_states*n_inf_states + 
                                                     n_exp_states*HH_5 + 
                                                     n_inf_states*HH_5 +
                                                     r_tau + r_beta,
                                                   data = df_control_measures_bias_vax_propH_effH)
summary(fit_hh_peak_size_bias_rel_int_vax_propH_effH)

fit_hh_epidemic_size_rel_int_vax_propH_effH <- lm(CumInfTot_diff_bias_rel ~ # n_exp_states*n_inf_states + 
                                                    n_exp_states*HH_5 + 
                                                    n_inf_states*HH_5 + 
                                                    r_tau + r_beta, 
                                                  data = df_control_measures_bias_vax_propH_effH)
summary(fit_hh_epidemic_size_rel_int_vax_propH_effH)
# Linear combination coefficient tests
l_lincomm_rel_int_vax_propH_effH <- vector(mode = "list", length = 3)
l_lincomm_rel_int_vax_propH_effH[[1]] <- summary(multcomp::glht(fit_hh_peak_time_bias_rel_int_vax_propH_effH, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))
l_lincomm_rel_int_vax_propH_effH[[2]] <- summary(multcomp::glht(fit_hh_peak_size_bias_rel_int_vax_propH_effH, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))
l_lincomm_rel_int_vax_propH_effH[[3]] <- summary(multcomp::glht(fit_hh_epidemic_size_rel_int_vax_propH_effH, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))

df_lincomm_rel_int_vax_propH_effH <- data.frame(`Epidemic outcome` = c("Peak time", "Peak size", "Epidemic size"),
                                                coef = as.numeric(c(l_lincomm_rel_int_vax_propH_effH[[1]]$test$coefficients, 
                                                                    l_lincomm_rel_int_vax_propH_effH[[2]]$test$coefficients,
                                                                    l_lincomm_rel_int_vax_propH_effH[[3]]$test$coefficients)), 
                                                pval = c(l_lincomm_rel_int_vax_propH_effH[[1]]$test$pvalues[1], 
                                                         l_lincomm_rel_int_vax_propH_effH[[2]]$test$pvalues[1],
                                                         l_lincomm_rel_int_vax_propH_effH[[3]]$test$pvalues[1]), 
                                                check.names = F)
df_lincomm_rel_int_vax_propH_effH
xtable::xtable(df_lincomm_rel_int_vax_propH_effH)

df_lincomm_rel_vax_propH_effH$Interaction <- "No"
df_lincomm_rel_int_vax_propH_effH$Interaction <- "Yes"

df_lincomm_rel_vax_propH_effH_comb <- bind_rows(df_lincomm_rel_vax_propH_effH, 
                                                df_lincomm_rel_int_vax_propH_effH) %>%
  group_by(`Epidemic outcome`) %>%
  arrange(desc(`Epidemic outcome`), Interaction)
df_lincomm_rel_vax_propH_effH_comb

v_pval_rel_vax_propH_effH <- df_lincomm_rel_vax_propH_effH_comb$pval
v_star_rel_vax_propH_effH <- cut(v_pval_rel_vax_propH_effH,
                                 breaks = c(-Inf, 0.01, 0.05, 0.1, Inf), 
                                 labels = c("$^{***}$", "$^{**}$", "$^{*}$", ""))

v_coef_rel_vax_propH_effH <- as.character(round(df_lincomm_rel_vax_propH_effH_comb$coef, 3))
v_coef_sg_rel_vax_propH_effH <- paste0(v_coef_rel_vax_propH_effH, v_star_rel_vax_propH_effH)

stargazer(fit_hh_peak_time_bias_rel_vax_propH_effH, fit_hh_peak_time_bias_rel_int_vax_propH_effH,
          fit_hh_peak_size_bias_rel_vax_propH_effH, fit_hh_peak_size_bias_rel_int_vax_propH_effH,
          fit_hh_epidemic_size_rel_vax_propH_effH, fit_hh_epidemic_size_rel_int_vax_propH_effH,
          type = "text", 
          title = "Metaregression estimates vaccination coverage = 90% & effectiveness = 90%",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          label = "tab:vax_metareg_90_90",
          order = v_stargazer_int_effects_order_main,
          covariate.labels = v_stargazer_int_effects_label_main,
          add.lines = list(c("HH3 + HH5", 
                             v_coef_sg_rel_vax_propH_effH)),
          notes = v_stargazer_notes_int_effects,
          notes.append = TRUE,
          notes.align = "l",
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE)

stargazer(fit_hh_peak_time_bias_rel_vax_propH_effH, fit_hh_peak_time_bias_rel_int_vax_propH_effH,
          fit_hh_peak_size_bias_rel_vax_propH_effH, fit_hh_peak_size_bias_rel_int_vax_propH_effH,
          fit_hh_epidemic_size_rel_vax_propH_effH, fit_hh_epidemic_size_rel_int_vax_propH_effH,
          type = "latex", 
          out = "output/TableA3_metaregression_bias_vax_propH_effH.tex", 
          title = "Metaregression estimates vaccination coverage = 90% & effectiveness = 90%",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          label = "tab:vax_metareg_90_90",
          order = v_stargazer_int_effects_order_main,
          covariate.labels = v_stargazer_int_effects_label_main,
          add.lines = list(c("HH3 + HH5", 
                             v_coef_sg_rel_vax_propH_effH)),
          notes = v_stargazer_notes_int_effects,
          notes.append = TRUE,
          notes.align = "l",
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE)

### Table A4: Meta regression on category-specific outcomes on Vaccine Measures Coverage = 90% and Effectiveness = 50% NPI 0% reduction ----
df_control_measures_bias_vax_propH_effL <- df_out_inf_all_control_summ_ref_bias %>% 
  filter(r_omega == 0.000 & 
           !((`Multicompartment structure` %in% c("E=1, I=1") & n_hhsize %in% c(1))) & 
           eff_vax %in% c(0.5) & vax_prop == 0.9, level_npi %in% c(1.0)) 

## Without interactions
fit_hh_peak_time_bias_rel_vax_propH_effL <- lm(max_Inftot_time_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                                 HH_5 + 
                                                 r_tau + r_beta, 
                                               data = df_control_measures_bias_vax_propH_effL %>% 
                                                 filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel_vax_propH_effL)

fit_hh_peak_size_bias_rel_vax_propH_effL <- lm(max_Inftot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                                 HH_5 + 
                                                 r_tau + r_beta,
                                               data = df_control_measures_bias_vax_propH_effL)
summary(fit_hh_peak_size_bias_rel_vax_propH_effL)

fit_hh_epidemic_size_rel_vax_propH_effL <- lm(CumInfTot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                                HH_5 + 
                                                r_tau + r_beta, 
                                              data = df_control_measures_bias_vax_propH_effL)
summary(fit_hh_epidemic_size_rel_vax_propH_effL)
# Linear combination coefficient tests
l_lincomm_rel_vax_propH_effL <- vector(mode = "list", length = 3)
l_lincomm_rel_vax_propH_effL[[1]] <- summary(multcomp::glht(fit_hh_peak_time_bias_rel_vax_propH_effL, linfct = c("(Intercept) + HH_5TRUE = 0")))
l_lincomm_rel_vax_propH_effL[[2]] <- summary(multcomp::glht(fit_hh_peak_size_bias_rel_vax_propH_effL, linfct = c("(Intercept) + HH_5TRUE = 0")))
l_lincomm_rel_vax_propH_effL[[3]] <- summary(multcomp::glht(fit_hh_epidemic_size_rel_vax_propH_effL, linfct = c("(Intercept) + HH_5TRUE = 0")))

df_lincomm_rel_vax_propH_effL <- data.frame(`Epidemic outcome` = c("Peak time", "Peak size", "Epidemic size"),
                                            coef = as.numeric(c(l_lincomm_rel_vax_propH_effL[[1]]$test$coefficients, 
                                                                l_lincomm_rel_vax_propH_effL[[2]]$test$coefficients,
                                                                l_lincomm_rel_vax_propH_effL[[3]]$test$coefficients)), 
                                            pval = c(l_lincomm_rel_vax_propH_effL[[1]]$test$pvalues[1], 
                                                     l_lincomm_rel_vax_propH_effL[[2]]$test$pvalues[1],
                                                     l_lincomm_rel_vax_propH_effL[[3]]$test$pvalues[1]), 
                                            check.names = F)
df_lincomm_rel_vax_propH_effL
xtable::xtable(df_lincomm_rel_vax_propH_effL)

## With interactions
fit_hh_peak_time_bias_rel_int_vax_propH_effL <- lm(max_Inftot_time_diff_bias_rel ~ # n_exp_states*n_inf_states + 
                                                     n_exp_states*HH_5 + 
                                                     n_inf_states*HH_5 + 
                                                     r_tau + r_beta, 
                                                   data = df_control_measures_bias_vax_propH_effL %>% 
                                                     filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel_int_vax_propH_effL)

fit_hh_peak_size_bias_rel_int_vax_propH_effL <- lm(max_Inftot_diff_bias_rel ~ # n_exp_states*n_inf_states + 
                                                     n_exp_states*HH_5 + 
                                                     n_inf_states*HH_5 +
                                                     r_tau + r_beta,
                                                   data = df_control_measures_bias_vax_propH_effL)
summary(fit_hh_peak_size_bias_rel_int_vax_propH_effL)

fit_hh_epidemic_size_rel_int_vax_propH_effL <- lm(CumInfTot_diff_bias_rel ~ # n_exp_states*n_inf_states + 
                                                    n_exp_states*HH_5 + 
                                                    n_inf_states*HH_5 + 
                                                    r_tau + r_beta, 
                                                  data = df_control_measures_bias_vax_propH_effL)
summary(fit_hh_epidemic_size_rel_int_vax_propH_effL)
# Linear combination coefficient tests
l_lincomm_rel_int_vax_propH_effL <- vector(mode = "list", length = 3)
l_lincomm_rel_int_vax_propH_effL[[1]] <- summary(multcomp::glht(fit_hh_peak_time_bias_rel_int_vax_propH_effL, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))
l_lincomm_rel_int_vax_propH_effL[[2]] <- summary(multcomp::glht(fit_hh_peak_size_bias_rel_int_vax_propH_effL, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))
l_lincomm_rel_int_vax_propH_effL[[3]] <- summary(multcomp::glht(fit_hh_epidemic_size_rel_int_vax_propH_effL, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))

df_lincomm_rel_int_vax_propH_effL <- data.frame(`Epidemic outcome` = c("Peak time", "Peak size", "Epidemic size"),
                                                coef = as.numeric(c(l_lincomm_rel_int_vax_propH_effL[[1]]$test$coefficients, 
                                                                    l_lincomm_rel_int_vax_propH_effL[[2]]$test$coefficients,
                                                                    l_lincomm_rel_int_vax_propH_effL[[3]]$test$coefficients)), 
                                                pval = c(l_lincomm_rel_int_vax_propH_effL[[1]]$test$pvalues[1], 
                                                         l_lincomm_rel_int_vax_propH_effL[[2]]$test$pvalues[1],
                                                         l_lincomm_rel_int_vax_propH_effL[[3]]$test$pvalues[1]), 
                                                check.names = F)
df_lincomm_rel_int_vax_propH_effL
xtable::xtable(df_lincomm_rel_int_vax_propH_effL)

df_lincomm_rel_vax_propH_effL$Interaction <- "No"
df_lincomm_rel_int_vax_propH_effL$Interaction <- "Yes"

df_lincomm_rel_vax_propH_effL_comb <- bind_rows(df_lincomm_rel_vax_propH_effL, 
                                                df_lincomm_rel_int_vax_propH_effL) %>%
  group_by(`Epidemic outcome`) %>%
  arrange(desc(`Epidemic outcome`), Interaction)
df_lincomm_rel_vax_propH_effL_comb

v_pval_rel_vax_propH_effL <- df_lincomm_rel_vax_propH_effL_comb$pval
v_star_rel_vax_propH_effL <- cut(v_pval_rel_vax_propH_effL,
                                 breaks = c(-Inf, 0.01, 0.05, 0.1, Inf), 
                                 labels = c("$^{***}$", "$^{**}$", "$^{*}$", ""))

v_coef_rel_vax_propH_effL <- as.character(round(df_lincomm_rel_vax_propH_effL_comb$coef, 3))
v_coef_sg_rel_vax_propH_effL <- paste0(v_coef_rel_vax_propH_effL, v_star_rel_vax_propH_effL)

stargazer(fit_hh_peak_time_bias_rel_vax_propH_effL, fit_hh_peak_time_bias_rel_int_vax_propH_effL,
          fit_hh_peak_size_bias_rel_vax_propH_effL, fit_hh_peak_size_bias_rel_int_vax_propH_effL,
          fit_hh_epidemic_size_rel_vax_propH_effL, fit_hh_epidemic_size_rel_int_vax_propH_effL,
          type = "text", 
          title = "Metaregression estimates vaccination coverage = 90% & effectiveness = 50%",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          label = "tab:vax_metareg_90_50",
          order = v_stargazer_int_effects_order_main,
          covariate.labels = v_stargazer_int_effects_label_main,
          add.lines = list(c("HH3 + HH5", 
                             v_coef_sg_rel_vax_propH_effL)),
          notes = v_stargazer_notes_int_effects,
          notes.append = TRUE,
          notes.align = "l",
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE)

stargazer(fit_hh_peak_time_bias_rel_vax_propH_effL, fit_hh_peak_time_bias_rel_int_vax_propH_effL,
          fit_hh_peak_size_bias_rel_vax_propH_effL, fit_hh_peak_size_bias_rel_int_vax_propH_effL,
          fit_hh_epidemic_size_rel_vax_propH_effL, fit_hh_epidemic_size_rel_int_vax_propH_effL,
          type = "latex", 
          out = "output/TableA4_metaregression_bias_vax_propH_effL.tex", 
          title = "Metaregression estimates vaccination coverage = 90% & effectiveness = 50%",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          label = "tab:vax_metareg_90_50",
          order = v_stargazer_int_effects_order_main,
          covariate.labels = v_stargazer_int_effects_label_main,
          add.lines = list(c("HH3 + HH5", 
                             v_coef_sg_rel_vax_propH_effL)),
          notes = v_stargazer_notes_int_effects,
          notes.append = TRUE,
          notes.align = "l",
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE)

### Table A5: Meta regression on category-specific outcomes on Vaccine Measures Coverage = 30% and Effectiveness = 50% NPI 0% reduction ----
df_control_measures_bias_vax_propL_effL <- df_out_inf_all_control_summ_ref_bias %>% 
  filter(r_omega == 0.000 & 
           !((`Multicompartment structure` %in% c("E=1, I=1") & n_hhsize %in% c(1))) & 
           eff_vax %in% c(0.5) & vax_prop == 0.3, level_npi %in% c(1.0)) 

## Without interactions
fit_hh_peak_time_bias_rel_vax_propL_effL <- lm(max_Inftot_time_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                                 HH_5 + 
                                                 r_tau + r_beta, 
                                               data = df_control_measures_bias_vax_propL_effL %>% 
                                                 filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel_vax_propL_effL)

fit_hh_peak_size_bias_rel_vax_propL_effL <- lm(max_Inftot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                                 HH_5 + 
                                                 r_tau + r_beta,
                                               data = df_control_measures_bias_vax_propL_effL)
summary(fit_hh_peak_size_bias_rel_vax_propL_effL)

fit_hh_epidemic_size_rel_vax_propL_effL <- lm(CumInfTot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                                HH_5 + 
                                                r_tau + r_beta, 
                                              data = df_control_measures_bias_vax_propL_effL)
summary(fit_hh_epidemic_size_rel_vax_propL_effL)
# Linear combination coefficient tests
l_lincomm_rel_vax_propL_effL <- vector(mode = "list", length = 3)
l_lincomm_rel_vax_propL_effL[[1]] <- summary(multcomp::glht(fit_hh_peak_time_bias_rel_vax_propL_effL, linfct = c("(Intercept) + HH_5TRUE = 0")))
l_lincomm_rel_vax_propL_effL[[2]] <- summary(multcomp::glht(fit_hh_peak_size_bias_rel_vax_propL_effL, linfct = c("(Intercept) + HH_5TRUE = 0")))
l_lincomm_rel_vax_propL_effL[[3]] <- summary(multcomp::glht(fit_hh_epidemic_size_rel_vax_propL_effL, linfct = c("(Intercept) + HH_5TRUE = 0")))

df_lincomm_rel_vax_propL_effL <- data.frame(`Epidemic outcome` = c("Peak time", "Peak size", "Epidemic size"),
                                            coef = as.numeric(c(l_lincomm_rel_vax_propL_effL[[1]]$test$coefficients, 
                                                                l_lincomm_rel_vax_propL_effL[[2]]$test$coefficients,
                                                                l_lincomm_rel_vax_propL_effL[[3]]$test$coefficients)), 
                                            pval = c(l_lincomm_rel_vax_propL_effL[[1]]$test$pvalues[1], 
                                                     l_lincomm_rel_vax_propL_effL[[2]]$test$pvalues[1],
                                                     l_lincomm_rel_vax_propL_effL[[3]]$test$pvalues[1]), 
                                            check.names = F)
df_lincomm_rel_vax_propL_effL
xtable::xtable(df_lincomm_rel_vax_propL_effL)

## With interactions
fit_hh_peak_time_bias_rel_int_vax_propL_effL <- lm(max_Inftot_time_diff_bias_rel ~ # n_exp_states*n_inf_states + 
                                                     n_exp_states*HH_5 + 
                                                     n_inf_states*HH_5 + 
                                                     r_tau + r_beta, 
                                                   data = df_control_measures_bias_vax_propL_effL %>% 
                                                     filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel_int_vax_propL_effL)

fit_hh_peak_size_bias_rel_int_vax_propL_effL <- lm(max_Inftot_diff_bias_rel ~ # n_exp_states*n_inf_states + 
                                                     n_exp_states*HH_5 + 
                                                     n_inf_states*HH_5 +
                                                     r_tau + r_beta,
                                                   data = df_control_measures_bias_vax_propL_effL)
summary(fit_hh_peak_size_bias_rel_int_vax_propL_effL)

fit_hh_epidemic_size_rel_int_vax_propL_effL <- lm(CumInfTot_diff_bias_rel ~ # n_exp_states*n_inf_states + 
                                                    n_exp_states*HH_5 + 
                                                    n_inf_states*HH_5 + 
                                                    r_tau + r_beta, 
                                                  data = df_control_measures_bias_vax_propL_effL)
summary(fit_hh_epidemic_size_rel_int_vax_propL_effL)
# Linear combination coefficient tests
l_lincomm_rel_int_vax_propL_effL <- vector(mode = "list", length = 3)
l_lincomm_rel_int_vax_propL_effL[[1]] <- summary(multcomp::glht(fit_hh_peak_time_bias_rel_int_vax_propL_effL, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))
l_lincomm_rel_int_vax_propL_effL[[2]] <- summary(multcomp::glht(fit_hh_peak_size_bias_rel_int_vax_propL_effL, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))
l_lincomm_rel_int_vax_propL_effL[[3]] <- summary(multcomp::glht(fit_hh_epidemic_size_rel_int_vax_propL_effL, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))

df_lincomm_rel_int_vax_propL_effL <- data.frame(`Epidemic outcome` = c("Peak time", "Peak size", "Epidemic size"),
                                                coef = as.numeric(c(l_lincomm_rel_int_vax_propL_effL[[1]]$test$coefficients, 
                                                                    l_lincomm_rel_int_vax_propL_effL[[2]]$test$coefficients,
                                                                    l_lincomm_rel_int_vax_propL_effL[[3]]$test$coefficients)), 
                                                pval = c(l_lincomm_rel_int_vax_propL_effL[[1]]$test$pvalues[1], 
                                                         l_lincomm_rel_int_vax_propL_effL[[2]]$test$pvalues[1],
                                                         l_lincomm_rel_int_vax_propL_effL[[3]]$test$pvalues[1]), 
                                                check.names = F)
df_lincomm_rel_int_vax_propL_effL
xtable::xtable(df_lincomm_rel_int_vax_propL_effL)

df_lincomm_rel_vax_propL_effL$Interaction <- "No"
df_lincomm_rel_int_vax_propL_effL$Interaction <- "Yes"

df_lincomm_rel_vax_propL_effL_comb <- bind_rows(df_lincomm_rel_vax_propL_effL, 
                                                df_lincomm_rel_int_vax_propL_effL) %>%
  group_by(`Epidemic outcome`) %>%
  arrange(desc(`Epidemic outcome`), Interaction)
df_lincomm_rel_vax_propL_effL_comb

v_pval_rel_vax_propL_effL <- df_lincomm_rel_vax_propL_effL_comb$pval
v_star_rel_vax_propL_effL <- cut(v_pval_rel_vax_propL_effL,
                                 breaks = c(-Inf, 0.01, 0.05, 0.1, Inf), 
                                 labels = c("$^{***}$", "$^{**}$", "$^{*}$", ""))

v_coef_rel_vax_propL_effL <- as.character(round(df_lincomm_rel_vax_propL_effL_comb$coef, 3))
v_coef_sg_rel_vax_propL_effL <- paste0(v_coef_rel_vax_propL_effL, v_star_rel_vax_propH_effL)

stargazer(fit_hh_peak_time_bias_rel_vax_propL_effL, fit_hh_peak_time_bias_rel_int_vax_propL_effL,
          fit_hh_peak_size_bias_rel_vax_propL_effL, fit_hh_peak_size_bias_rel_int_vax_propL_effL,
          fit_hh_epidemic_size_rel_vax_propL_effL, fit_hh_epidemic_size_rel_int_vax_propL_effL,
          type = "text", 
          title = "Metaregression estimates vaccination coverage = 30% & effectiveness = 50%",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          label = "tab:vax_metareg_30_50",
          order = v_stargazer_int_effects_order_main,
          covariate.labels = v_stargazer_int_effects_label_main,
          add.lines = list(c("HH3 + HH5", 
                             v_coef_sg_rel_vax_propL_effL)),
          notes = v_stargazer_notes_int_effects,
          notes.append = TRUE,
          notes.align = "l",
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE)

stargazer(fit_hh_peak_time_bias_rel_vax_propL_effL, fit_hh_peak_time_bias_rel_int_vax_propL_effL,
          fit_hh_peak_size_bias_rel_vax_propL_effL, fit_hh_peak_size_bias_rel_int_vax_propL_effL,
          fit_hh_epidemic_size_rel_vax_propL_effL, fit_hh_epidemic_size_rel_int_vax_propL_effL,
          type = "latex", 
          out = "output/TableA5_metaregression_bias_vax_propL_effL.tex", 
          title = "Metaregression estimates vaccination coverage = 30% & effectiveness = 50%",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          label = "tab:vax_metareg_30_50",
          order = v_stargazer_int_effects_order_main,
          covariate.labels = v_stargazer_int_effects_label_main,
          add.lines = list(c("HH3 + HH5", 
                             v_coef_sg_rel_vax_propL_effL)),
          notes = v_stargazer_notes_int_effects,
          notes.append = TRUE,
          notes.align = "l",
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE)

### Table A6: Meta regression on category-specific outcomes on Control Measures NPI 60% reduction, OMEGA = 0.01 ----
df_out_inf_all_control_summ_ref_bias$HH_5 <- df_out_inf_all_control_summ_ref_bias$n_hhsize == 5
df_out_inf_all_control_summ_ref_bias_metareg_W01 <- df_out_inf_all_control_summ_ref_bias %>% 
  filter(r_omega == 0.01 & 
           !((`Multicompartment structure` %in% c("E=1, I=1") & n_hhsize %in% c(1))) & 
           eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi %in% c(0.4)) 

## Without interactions
fit_hh_peak_time_bias_rel_W01 <- lm(max_Inftot_time_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                  HH_5 +
                                  r_tau + r_beta, 
                                data = df_out_inf_all_control_summ_ref_bias_metareg_W01 %>% 
                                  filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel_W01)

fit_hh_peak_size_bias_rel_W01 <- lm(max_Inftot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                  HH_5 +
                                  r_tau + r_beta,
                                data = df_out_inf_all_control_summ_ref_bias_metareg_W01)
summary(fit_hh_peak_size_bias_rel_W01)

fit_hh_epidemic_size_rel_W01 <- lm(CumInfTot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                 HH_5 +
                                 r_tau + r_beta, 
                               data = df_out_inf_all_control_summ_ref_bias_metareg_W01)
summary(fit_hh_epidemic_size_rel_W01)

# Linear combination coefficient tests
l_lincomm_rel_W01 <- vector(mode = "list", length = 3)
l_lincomm_rel_W01[[1]] <- summary(multcomp::glht(fit_hh_peak_time_bias_rel_W01, linfct = c("(Intercept) + HH_5TRUE = 0")))
l_lincomm_rel_W01[[2]] <- summary(multcomp::glht(fit_hh_peak_size_bias_rel_W01, linfct = c("(Intercept) + HH_5TRUE = 0")))
l_lincomm_rel_W01[[3]] <- summary(multcomp::glht(fit_hh_epidemic_size_rel_W01, linfct = c("(Intercept) + HH_5TRUE = 0")))

df_lincomm_rel_W01 <- data.frame(`Epidemic outcome` = c("Peak time", "Peak size", "Epidemic size"),
                             coef = as.numeric(c(l_lincomm_rel_W01[[1]]$test$coefficients, 
                                                 l_lincomm_rel_W01[[2]]$test$coefficients,
                                                 l_lincomm_rel_W01[[3]]$test$coefficients)), 
                             pval = c(l_lincomm_rel_W01[[1]]$test$pvalues[1], 
                                      l_lincomm_rel_W01[[2]]$test$pvalues[1],
                                      l_lincomm_rel_W01[[3]]$test$pvalues[1]), 
                             check.names = F)
df_lincomm_rel_W01
xtable::xtable(df_lincomm_rel_W01)

## With interactions
fit_hh_peak_time_bias_rel_int_W01 <- lm(max_Inftot_time_diff_bias_rel ~ #n_exp_states*n_inf_states + 
                                      n_exp_states*HH_5 + 
                                      n_inf_states*HH_5 +
                                      r_tau + r_beta, 
                                    data = df_out_inf_all_control_summ_ref_bias_metareg_W01 %>% 
                                      filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel_int_W01)

fit_hh_peak_size_bias_rel_int_W01 <- lm(max_Inftot_diff_bias_rel ~ #n_exp_states*n_inf_states + 
                                      n_exp_states*HH_5 + 
                                      n_inf_states*HH_5 +
                                      r_tau + r_beta,
                                    data = df_out_inf_all_control_summ_ref_bias_metareg_W01)
summary(fit_hh_peak_size_bias_rel_int_W01)

fit_hh_epidemic_size_rel_int_W01 <- lm(CumInfTot_diff_bias_rel ~ #n_exp_states*n_inf_states + 
                                     n_exp_states*HH_5 + 
                                     n_inf_states*HH_5 +
                                     r_tau + r_beta, 
                                   data = df_out_inf_all_control_summ_ref_bias_metareg_W01)
summary(fit_hh_epidemic_size_rel_int_W01)

# Linear combination coefficient tests
l_lincomm_rel_int_W01 <- vector(mode = "list", length = 3)
l_lincomm_rel_int_W01[[1]] <- summary(multcomp::glht(fit_hh_peak_time_bias_rel_int_W01, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))
l_lincomm_rel_int_W01[[2]] <- summary(multcomp::glht(fit_hh_peak_size_bias_rel_int_W01, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))
l_lincomm_rel_int_W01[[3]] <- summary(multcomp::glht(fit_hh_epidemic_size_rel_int_W01, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))

df_lincomm_rel_int_W01 <- data.frame(`Epidemic outcome` = c("Peak time", "Peak size", "Epidemic size"),
                                 coef = as.numeric(c(l_lincomm_rel_int_W01[[1]]$test$coefficients, 
                                                     l_lincomm_rel_int_W01[[2]]$test$coefficients,
                                                     l_lincomm_rel_int_W01[[3]]$test$coefficients)), 
                                 pval = c(l_lincomm_rel_int_W01[[1]]$test$pvalues[1], 
                                          l_lincomm_rel_int_W01[[2]]$test$pvalues[1],
                                          l_lincomm_rel_int_W01[[3]]$test$pvalues[1]), 
                                 check.names = F)
df_lincomm_rel_int_W01
xtable::xtable(df_lincomm_rel_int_W01)

df_lincomm_rel_W01$Interaction <- "No"
df_lincomm_rel_int_W01$Interaction <- "Yes"

df_lincomm_rel_W01_comb <- bind_rows(df_lincomm_rel_W01, 
                                                df_lincomm_rel_int_W01) %>%
  group_by(`Epidemic outcome`) %>%
  arrange(desc(`Epidemic outcome`), Interaction)
df_lincomm_rel_W01_comb

v_pval_rel_W01 <- df_lincomm_rel_W01_comb$pval
v_star_rel_W01 <- cut(v_pval_rel_W01,
                                 breaks = c(-Inf, 0.01, 0.05, 0.1, Inf), 
                                 labels = c("$^{***}$", "$^{**}$", "$^{*}$", ""))

v_coef_rel_W01 <- as.character(round(df_lincomm_rel_W01_comb$coef, 3))
v_coef_sg_rel_W01 <- paste0(v_coef_rel_W01, v_star_rel_W01)

stargazer(fit_hh_peak_time_bias_rel_W01, fit_hh_peak_time_bias_rel_int_W01,
          fit_hh_peak_size_bias_rel_W01, fit_hh_peak_size_bias_rel_int_W01,
          fit_hh_epidemic_size_rel_W01, fit_hh_epidemic_size_rel_int_W01,
          type = "text", 
          title = "Metaregression estimates on relative bias of treatment effects, NPI = 60% and $\\omega$ = 0.01",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          label = "tab:vax_metareg_NPI60_W01",
          order = v_stargazer_int_effects_order_main,
          covariate.labels = v_stargazer_int_effects_label_main,
          add.lines = list(c("HH3 + HH5", 
                             v_coef_sg_rel_W01)),
          notes = v_stargazer_notes_int_effects,
          notes.append = TRUE,
          notes.align = "l",
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE)

stargazer(fit_hh_peak_time_bias_rel_W01, fit_hh_peak_time_bias_rel_int_W01,
          fit_hh_peak_size_bias_rel_W01, fit_hh_peak_size_bias_rel_int_W01,
          fit_hh_epidemic_size_rel_W01, fit_hh_epidemic_size_rel_int_W01,
          type = "latex", 
          out = "output/TableA6_metaregression_bias_W01.tex", 
          title = "Metaregression estimates on relative bias of treatment effects, NPI = 60% and $\\omega$ = 0.01",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          label = "tab:vax_metareg_NPI60_W01",
          order = v_stargazer_int_effects_order_main,
          covariate.labels = v_stargazer_int_effects_label_main,
          add.lines = list(c("HH3 + HH5", 
                             v_coef_sg_rel_W01)),
          notes = v_stargazer_notes_int_effects,
          notes.append = TRUE,
          notes.align = "l",
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE)

### Table A7: Meta regression on category-specific outcomes on Control Measures NPI 60% reduction, OMEGA = 0.02 ----
df_out_inf_all_control_summ_ref_bias_metareg_W02 <- df_out_inf_all_control_summ_ref_bias %>% 
  filter(r_omega == 0.02 & 
           !((`Multicompartment structure` %in% c("E=1, I=1") & n_hhsize %in% c(1))) & 
           eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi %in% c(0.4)) 

## Without interactions
fit_hh_peak_time_bias_rel_W02 <- lm(max_Inftot_time_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                      HH_5 +
                                      r_tau + r_beta, 
                                    data = df_out_inf_all_control_summ_ref_bias_metareg_W02 %>% 
                                      filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel_W02)

fit_hh_peak_size_bias_rel_W02 <- lm(max_Inftot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                      HH_5 +
                                      r_tau + r_beta,
                                    data = df_out_inf_all_control_summ_ref_bias_metareg_W02)
summary(fit_hh_peak_size_bias_rel_W02)

fit_hh_epidemic_size_rel_W02 <- lm(CumInfTot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                     HH_5 +
                                     r_tau + r_beta, 
                                   data = df_out_inf_all_control_summ_ref_bias_metareg_W02)
summary(fit_hh_epidemic_size_rel_W02)

# Linear combination coefficient tests
l_lincomm_rel_W02 <- vector(mode = "list", length = 3)
l_lincomm_rel_W02[[1]] <- summary(multcomp::glht(fit_hh_peak_time_bias_rel_W02, linfct = c("(Intercept) + HH_5TRUE = 0")))
l_lincomm_rel_W02[[2]] <- summary(multcomp::glht(fit_hh_peak_size_bias_rel_W02, linfct = c("(Intercept) + HH_5TRUE = 0")))
l_lincomm_rel_W02[[3]] <- summary(multcomp::glht(fit_hh_epidemic_size_rel_W02, linfct = c("(Intercept) + HH_5TRUE = 0")))

df_lincomm_rel_W02 <- data.frame(`Epidemic outcome` = c("Peak time", "Peak size", "Epidemic size"),
                                 coef = as.numeric(c(l_lincomm_rel_W02[[1]]$test$coefficients, 
                                                     l_lincomm_rel_W02[[2]]$test$coefficients,
                                                     l_lincomm_rel_W02[[3]]$test$coefficients)), 
                                 pval = c(l_lincomm_rel_W02[[1]]$test$pvalues[1], 
                                          l_lincomm_rel_W02[[2]]$test$pvalues[1],
                                          l_lincomm_rel_W02[[3]]$test$pvalues[1]), 
                                 check.names = F)
df_lincomm_rel_W02
xtable::xtable(df_lincomm_rel_W02)

## With interactions
fit_hh_peak_time_bias_rel_int_W02 <- lm(max_Inftot_time_diff_bias_rel ~ #n_exp_states*n_inf_states + 
                                          n_exp_states*HH_5 + 
                                          n_inf_states*HH_5 +
                                          r_tau + r_beta, 
                                        data = df_out_inf_all_control_summ_ref_bias_metareg_W02 %>% 
                                          filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel_int_W02)

fit_hh_peak_size_bias_rel_int_W02 <- lm(max_Inftot_diff_bias_rel ~ #n_exp_states*n_inf_states + 
                                          n_exp_states*HH_5 + 
                                          n_inf_states*HH_5 +
                                          r_tau + r_beta,
                                        data = df_out_inf_all_control_summ_ref_bias_metareg_W02)
summary(fit_hh_peak_size_bias_rel_int_W02)

fit_hh_epidemic_size_rel_int_W02 <- lm(CumInfTot_diff_bias_rel ~ #n_exp_states*n_inf_states + 
                                         n_exp_states*HH_5 + 
                                         n_inf_states*HH_5 +
                                         r_tau + r_beta, 
                                       data = df_out_inf_all_control_summ_ref_bias_metareg_W02)
summary(fit_hh_epidemic_size_rel_int_W02)

# Linear combination coefficient tests
l_lincomm_rel_int_W02 <- vector(mode = "list", length = 3)
l_lincomm_rel_int_W02[[1]] <- summary(multcomp::glht(fit_hh_peak_time_bias_rel_int_W02, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))
l_lincomm_rel_int_W02[[2]] <- summary(multcomp::glht(fit_hh_peak_size_bias_rel_int_W02, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))
l_lincomm_rel_int_W02[[3]] <- summary(multcomp::glht(fit_hh_epidemic_size_rel_int_W02, linfct = c("(Intercept) + HH_5TRUE + n_exp_states:HH_5TRUE + HH_5TRUE:n_inf_states = 0")))

df_lincomm_rel_int_W02 <- data.frame(`Epidemic outcome` = c("Peak time", "Peak size", "Epidemic size"),
                                     coef = as.numeric(c(l_lincomm_rel_int_W02[[1]]$test$coefficients, 
                                                         l_lincomm_rel_int_W02[[2]]$test$coefficients,
                                                         l_lincomm_rel_int_W02[[3]]$test$coefficients)), 
                                     pval = c(l_lincomm_rel_int_W02[[1]]$test$pvalues[1], 
                                              l_lincomm_rel_int_W02[[2]]$test$pvalues[1],
                                              l_lincomm_rel_int_W02[[3]]$test$pvalues[1]), 
                                     check.names = F)
df_lincomm_rel_int_W02
xtable::xtable(df_lincomm_rel_int_W02)

df_lincomm_rel_W02$Interaction <- "No"
df_lincomm_rel_int_W02$Interaction <- "Yes"

df_lincomm_rel_W02_comb <- bind_rows(df_lincomm_rel_W02, 
                                     df_lincomm_rel_int_W02) %>%
  group_by(`Epidemic outcome`) %>%
  arrange(desc(`Epidemic outcome`), Interaction)
df_lincomm_rel_W02_comb

v_pval_rel_W02 <- df_lincomm_rel_W02_comb$pval
v_star_rel_W02 <- cut(v_pval_rel_W02,
                      breaks = c(-Inf, 0.01, 0.05, 0.1, Inf), 
                      labels = c("$^{***}$", "$^{**}$", "$^{*}$", ""))

v_coef_rel_W02 <- as.character(round(df_lincomm_rel_W02_comb$coef, 3))
v_coef_sg_rel_W02 <- paste0(v_coef_rel_W02, v_star_rel_W02)

stargazer(fit_hh_peak_time_bias_rel_W02, fit_hh_peak_time_bias_rel_int_W02,
          fit_hh_peak_size_bias_rel_W02, fit_hh_peak_size_bias_rel_int_W02,
          fit_hh_epidemic_size_rel_W02, fit_hh_epidemic_size_rel_int_W02,
          type = "text", 
          title = "Metaregression estimates on relative bias of treatment effects, NPI = 60% and $\\omega$ = 0.02",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          label = "tab:vax_metareg_NPI60_W02",
          order = v_stargazer_int_effects_order_main,
          covariate.labels = v_stargazer_int_effects_label_main,
          add.lines = list(c("HH3 + HH5", 
                             v_coef_sg_rel_W02)),
          notes = v_stargazer_notes_int_effects,
          notes.append = TRUE,
          notes.align = "l",
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE)

stargazer(fit_hh_peak_time_bias_rel_W02, fit_hh_peak_time_bias_rel_int_W02,
          fit_hh_peak_size_bias_rel_W02, fit_hh_peak_size_bias_rel_int_W02,
          fit_hh_epidemic_size_rel_W02, fit_hh_epidemic_size_rel_int_W02,
          type = "latex", 
          out = "output/TableA7_metaregression_bias_W02.tex", 
          title = "Metaregression estimates on relative bias of treatment effects, NPI = 60% and $\\omega$ = 0.02",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          label = "tab:vax_metareg_NPI60_W02",
          order = v_stargazer_int_effects_order_main,
          covariate.labels = v_stargazer_int_effects_label_main,
          add.lines = list(c("HH3 + HH5", 
                             v_coef_sg_rel_W02)),
          notes = v_stargazer_notes_int_effects,
          notes.append = TRUE,
          notes.align = "l",
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE)

# Save workspace of all work up to this date ----
save.image(file = "output/WorkSpace - HHMCSEIRV - 2022-06-23.RData")
# load(file = "output/WorkSpace - HHMCSEIRV - 2022-06-23.RData")
