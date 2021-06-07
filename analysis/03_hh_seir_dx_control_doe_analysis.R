library(dplyr)
library(ggplot2)
source("R/02_decision_model_functions.R")


GLOBAL_PARAM_FILE <- "data/df_doe_mc_seirv_control.RDS"
df_doe_mc_seirv_runs <- readRDS(GLOBAL_PARAM_FILE)
v_pid <- unique(df_doe_mc_seirv_runs$pid)
n_doe <- length(v_pid)

list_l_out <- list()
v_failed_pid <- c()
for(i in 1:n_doe){ # n_pid <- v_pid[1]
  tryCatch({
    load(file = paste0("output/output_doe_mc_seirv_", v_pid[i],".RData")) 
    list_l_out[[i]] <- l_out
    names(list_l_out)[i] <- v_pid[i]
  }, error = function(cond) {
    print(paste(Sys.time(),": +++ ERROR in pid", v_pid[i], cond, "\n"))
    v_failed_pid <<- c(v_failed_pid, v_pid[i])
  })
  if(i/(n_doe/100) == round(i/(n_doe/100),0)) {
    cat('\r', paste(i/n_doe * 100, "% done", sep = " "))
  }
}

v_good_pid <- as.numeric(names(list_l_out))[!is.na(as.numeric(names(list_l_out)))]
# v_failed_pid <- v_pid[is.na(as.numeric(names(list_l_out)))]

save(list_l_out, v_good_pid,
     file = "output/output_doe_mc_seirv_all_control.RData")
gc()
#### Analyze output from DOE Natural History ####
load(file = "output/output_doe_mc_seirv_all_control.RData")

## Obtain runs for which DOE failed
df_doe_mc_seirv_runs_failed <- df_doe_mc_seirv_runs %>%
  filter(pid %in% v_failed_pid)

## compute epidemic outputs for each of the good pids
## Put all these outputs in a long data.frame
l_out_test <- list_l_out[[1]]
show_MC_SEIRV_model_results(l_out_test)
df_inf_test <- calc_inf_totals(l_out_test)

df_out_inf_all <- c()
for(n_pid in v_good_pid){ # n_pid <- v_good_pid[1]
  l_out_temp <- list_l_out[[as.character(n_pid)]]
  df_inf_temp <- data.frame(df_doe_mc_seirv_runs %>%
                              filter(pid %in% n_pid),
                            calc_inf_nodx(l_out_temp),
                            Inftot = calc_inf_totals(l_out_temp)$Inftot)
  df_out_inf_all <- bind_rows(df_out_inf_all, 
                              df_inf_temp)
  
}
save(df_out_inf_all, 
     file = "output/df_output_doe_mc_seirv_all_control.RData")

#### Analyze epidemic outputs ####
load(file = "output/df_output_doe_mc_seirv_all_control.RData")
df_out_inf_all$Esize <- paste0("# of E compartments = ", df_out_inf_all$n_exp_states)
df_out_inf_all$Isize <- paste0("# of I compartments = ", df_out_inf_all$n_inf_states)
df_out_inf_all$`Household size` <- ordered(df_out_inf_all$n_hhsize, unique(df_out_inf_all$n_hhsize))
df_out_inf_all$n_hhsize <- ordered(df_out_inf_all$n_hhsize)
df_out_inf_all$PropVax <- paste0("Proportion vaccinated = ", scales::percent(df_out_inf_all$vax_prop))
df_out_inf_all$NPIeff <- paste0("NPI effectiveness = ", scales::percent(1-df_out_inf_all$level_npi))

df_out_inf_all_summ <- df_out_inf_all %>%
  group_by(pid) %>%
  mutate(# Find the t at which I(t) is at its max
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

hist(df_out_inf_all_summ$max_Inftot)
hist(df_out_inf_all_summ$max_Inftot_time)
hist(df_out_inf_all_summ$IDX100_time)
hist(df_out_inf_all_summ$IDX500_time)
hist(df_out_inf_all_summ$p05_Inftot_time, breaks = 15)
hist(df_out_inf_all_summ$p10_Inftot_time, breaks = 15)
hist(df_out_inf_all_summ$p25_Inftot_time, breaks = 15)
hist(df_out_inf_all_summ$p50_Inftot_time, breaks = 15)

df_out_inf_all %>% 
  filter(r_beta == 0.25 & r_tau == 0.50 & r_omega == 0 & Inftot >=0)

df_out_inf_all %>% 
  filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.020 & time <= 70 & 
           n_hhsize == 5,
           n_exp_states == 3 & n_inf_states == 3,
         vax_prop == 0.9) %>%
  View()

ggplot(df_out_inf_all %>% 
         filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.020 & time <= 70 & level_npi == 1), 
       aes(x = time, y = Inftot, color = n_hhsize)) +
  geom_line(size = 1.3) +
  facet_grid(Esize ~ Isize) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom",
        strip.background = element_rect(colour="white", fill="lightgray"),
        # legend.position=c(.88,.3),
        legend.key = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin=margin(-10,-10,-10,-10))


ggplot(df_out_inf_all %>% 
         filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.000 & time <= 70 & 
                  n_exp_states == 3 & n_inf_states == 3 & n_hhsize %in% c(1, 3, 5) & 
                  eff_vax %in% c(0.9)), 
       aes(x = time, y = Inftot, color = `Household size`)) + # linetype = as.factor(eff_vax))
  geom_line(size = 1.1) +
  facet_grid(NPIeff ~ PropVax) +
  # scale_color_viridis_d(option = "C", direction = -1) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom",
        strip.background = element_rect(colour="white", fill="lightgray"),
        # legend.position=c(.88,.3),
        legend.key = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin=margin(-10,-10,-10,-10))

ggplot(df_out_inf_all_summ %>% 
         filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.000 & time <= 70 & 
                  n_exp_states == 3 & n_inf_states == 3 & n_hhsize %in% c(1, 3, 5) & 
                  eff_vax %in% c(0.9)), 
       aes(x = `Household size`, y = CumInfTot, fill = `Household size`)) + # linetype = as.factor(eff_vax))
  geom_col(color = NA) +
  facet_grid(NPIeff ~ PropVax) +
  scale_y_continuous(labels = function(x)round(x/10e6, digits = 2)) +
  # scale_fill_viridis_d(option = "C", direction = -1) +
  # scale_color_jcolors(palette = "rainbow") +
  # scale_fill_brewer(palette = "Spectral") +
  # xlab("Household size") +
  ylab("Cumulative infections (millions)") +
  theme_bw(base_size = 16) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        legend.position = c(""),
        # legend.position = c(.88,.8),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.key = element_blank())

ggplot(df_out_inf_all_summ %>% 
         filter(r_beta == 0.25 & r_tau == 0.50 & r_omega == 0.000 & time <= 70 & 
                  n_exp_states == 3 & n_inf_states == 3 & n_hhsize %in% c(1, 3, 5) & 
                  eff_vax %in% c(0.9)), 
       aes(x = `Household size`, y = max_Inftot, fill = `Household size`)) + # linetype = as.factor(eff_vax))
  geom_col(color = NA) +
  facet_grid(NPIeff ~ PropVax) +
  scale_y_continuous(labels = function(x) scales::percent(x/10e6, accuracy = 1.0)) +
  # scale_fill_viridis_d(option = "D", direction = -1) +
  scale_fill_grey() +
  # scale_color_jcolors(palette = "rainbow") +
  # scale_fill_brewer(palette = "Spectral") +
  # xlab("Household size") +
  ylab("Magnitude of epidemic peak (% of total population)") +
  theme_bw(base_size = 16) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        legend.position = c(""),
        # legend.position = c(.88,.8),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.key = element_blank())
ggsave("figs/SMDM_household_communiy_MC_SEIR_NPI_VAX.pdf", width = 11, height = 8)
ggsave("figs/SMDM_household_communiy_MC_SEIR_NPI_VAX.png", width = 11, height = 8)
ggsave("figs/SMDM_household_communiy_MC_SEIR_NPI_VAX.jpeg", width = 11, height = 8)

ggplot(df_out_inf_all %>% 
         filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.020 & time <= 70 & 
                  n_exp_states == 3 & n_inf_states == 3 &
                  eff_vax == 0.5), 
       aes(x = time, y = as.factor(eff_vax), fill = Inftot)) +
  geom_tile() +
  facet_grid((1-level_npi) + n_hhsize ~ vax_prop) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")

# List of to-dos:
# - add computation time
# - Add Incident infections equation and capture that on the output
# HH vs non-HH transmission.

df_params_naming <- expand.grid(r_beta  = unique(df_out_inf_all$r_beta),
                                r_tau   = unique(df_out_inf_all$r_tau),
                                r_omega = unique(df_out_inf_all$r_omega))

for(i in 1:nrow(df_params_naming)){
  ggplot(df_out_inf_all %>% filter(r_beta == df_params_naming[i, ]$r_beta & 
                                     r_tau == df_params_naming[i, ]$r_tau & 
                                     r_omega ==df_params_naming[i, ]$r_omega & 
                                     time <= 70 & n_hhsize < 7), 
         aes(x = time, y = Inftot, color = n_hhsize)) +
    geom_line(size = 1.3) +
    facet_grid(n_exp_states ~ n_inf_states) +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom")
  ggsave(filename = paste0("figs/03_hh_seir_dx_nathist_doe_plot_", 
                           "o", df_params_naming[i, ]$r_omega, "_",
                           "t", df_params_naming[i, ]$r_tau, "_",
                           "b", df_params_naming[i, ]$r_beta, ".pdf"), 
         width = 8, height = 6)
}
# Ranges of times at which I starts to rise: compute doubling times
