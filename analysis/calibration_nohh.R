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


# Generate projections with calibrated betas in the absence of HH structure ----
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
for(i in 1:n_proj_nhh){ # i <- 1
  temp <- l_out_projection[[i]]
  v_names_exp    <- paste("E", letters[seq(1, temp$n_exp_states[1])], sep = "")
  v_names_inf    <- paste("I", letters[seq(1, temp$n_inf_states[1])], sep = "")
  v_names_inf_dx <- paste("IDX", letters[seq(1, temp$n_inf_states[1])], sep = "")
  
  df_temp <- data.frame(pid = temp$pid, 
                        level_npi = temp$level_npi,
                        eff_vax = temp$eff_vax,
                        vax_prop = temp$vax_prop,
                        time = temp$time,
                        Exptot_nhh = rowSums(temp[, v_names_exp,drop=FALSE]),
                        InfNoDX_nhh = rowSums(temp[, v_names_inf,drop=FALSE]),
                        Inftot_nhh = rowSums(temp[, c(v_names_inf,
                                                    v_names_inf_dx), drop=FALSE]),
                        check.names = FALSE)
  df_out_inf_all_nhh <- bind_rows(df_out_inf_all_nhh,
                                  df_temp)
}

df_out_inf_all_nhh_noint <- df_out_inf_all_nhh %>% 
  filter(vax_prop == 0,
         level_npi == 1)
# Rename variables for pretty plotting format
df_out_inf_all_nhh_noint$Esize <- paste0("# of E compartments = ", df_out_inf_all_nhh_noint$n_exp_states)
df_out_inf_all_nhh_noint$Isize <- paste0("# of I compartments = ", df_out_inf_all_nhh_noint$n_inf_states)


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
    IDX500_nhh_time = max(time[which((Inftot_nhh-InfNoDX_nhh) <= 500 & time < max_Inftot_nhh_time)]),
    IDX100_nhh_time = max(time[which((Inftot_nhh-InfNoDX_nhh) <= 100 & time < max_Inftot_nhh_time)]),
    CumInfTot_nhh = sum(Inftot_nhh)) %>%
  slice_head() %>%
  ungroup()

# Control measures outputs ----
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
  right_join(df_out_inf_all_summ %>% filter(n_hhsize > 1, 
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
  select(pid)

df_fig_nohh_vs_hh <-  bind_rows(df_out_inf_all_nhh_noint %>% 
                                  mutate(Structure = "No HH",
                                         Exptot = Exptot_nhh ,
                                         InfNoDX = InfNoDX_nhh ,
                                         Inftot = Inftot_nhh ), 
                                df_out_inf_all %>% mutate(Structure = "HH")) %>% 
  filter(pid %in% as.matrix(v_fig_pid))



gg_epidemic_curve_nohh_vs_hh <- ggplot(df_fig_nohh_vs_hh, 
                                                aes(x = time, y = Inftot/10e6, color = Structure)) + # 
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
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1.0), limits = c(0, 0.06)) +
  # scale_color_grey(start = 0.2, end = 0.6) +
  scale_color_manual(values = c("1" = "blue", "3" = "red")) +
  xlab("Time") +
  ylab("Infected population (% of total population)") +
  guides(color = guide_legend(nrow = 2), 
         linetype = guide_legend(nrow = 1)) +
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

# Metaregression ----
## Without interactions
fit_hh_peak_time_bias_rel_nhh <- lm(max_Inftot_time_diff_perc ~ n_exp_states + 
                                      n_inf_states + 
                                      n_hhsize + 
                                      r_tau + r_beta, 
                                    data = df_out_inf_all_biased_noint_summ %>% 
                                      filter(max_Inftot_time_diff_perc != Inf))
summary(fit_hh_peak_time_bias_rel_nhh)

fit_hh_peak_size_bias_rel_nhh <- lm(max_Inftot_diff_perc ~ n_exp_states + 
                                      n_inf_states + 
                                      n_hhsize + 
                                      r_tau + r_beta,
                                    data = df_out_inf_all_biased_noint_summ)
summary(fit_hh_peak_size_bias_rel_nhh)

fit_hh_epidemic_size_rel_nhh <- lm(CumInfTot_diff_perc ~ n_exp_states + 
                                     n_inf_states + 
                                     n_hhsize +
                                     r_tau + r_beta, 
                                   data = df_out_inf_all_biased_noint_summ)
summary(fit_hh_epidemic_size_rel_nhh)

## With interactions
fit_hh_peak_time_bias_rel_nhh <- lm(max_Inftot_time_diff_perc ~ n_exp_states*n_inf_states + 
                                      n_exp_states*n_hhsize + 
                                      n_inf_states*n_hhsize + 
                                      r_tau + r_beta, 
                                    data = df_out_inf_all_biased_noint_summ %>% 
                                      filter(max_Inftot_time_diff_perc != Inf))
summary(fit_hh_peak_time_bias_rel_nhh)

fit_hh_peak_size_bias_rel_nhh <- lm(max_Inftot_diff_perc ~ n_exp_states*n_inf_states + 
                                      n_exp_states*n_hhsize + 
                                      n_inf_states*n_hhsize +
                                      r_tau + r_beta,
                                    data = df_out_inf_all_biased_noint_summ)
summary(fit_hh_peak_size_bias_rel_nhh)

fit_hh_epidemic_size_rel_nhh <- lm(CumInfTot_diff_perc ~ n_exp_states*n_inf_states + 
                                     n_exp_states*n_hhsize + 
                                     n_inf_states*n_hhsize + 
                                     r_tau + r_beta, 
                                   data = df_out_inf_all_biased_noint_summ)
summary(fit_hh_epidemic_size_rel_nhh)

#* 1. summarize the biased projections
#* 2. Join in the true summarized projections with HH
#* 3. Compute biases 
#* 4. Visualize biases
#* 5. Metaregress biases
#* 6. Do similar things for intervention effects