rm(list = ls())
library(dplyr)
library(ggplot2)
library(patchwork)
library(stargazer)
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
# df_out_inf_all$n_hhsize <- ordered(df_out_inf_all$n_hhsize)
df_out_inf_all$`Vaccine effectiveness` <- scales::percent(df_out_inf_all$eff_vax)
df_out_inf_all$PropVax <- paste0("Proportion vaccinated = ", scales::percent(df_out_inf_all$vax_prop))
df_out_inf_all$NPIeff <- paste0("NPI effectiveness = ", scales::percent(1-df_out_inf_all$level_npi))
df_out_inf_all$`Multicompartment structure` <- paste0("E=", 
                                                      df_out_inf_all$n_exp_states, 
                                                      ", I=", 
                                                      df_out_inf_all$n_inf_states)

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

# Meta regression
fit_hh_peak_date <- lm(log(max_Inftot_time) ~ n_exp_states + n_inf_states + 
                         n_hhsize + r_tau + r_beta + r_omega, 
                       data = df_out_inf_all_summ)
summary(fit_hh_peak_date)
fit_hh_peak <- lm(log(max_Inftot) ~ n_exp_states + n_inf_states + 
                    n_hhsize + r_tau + r_beta + r_omega,
                  data = df_out_inf_all_summ)
summary(fit_hh_peak)
fit_hh_size <- lm(log(CumInfTot) ~ n_exp_states + n_inf_states +
                    n_hhsize + r_tau + r_beta + r_omega, 
                  data = df_out_inf_all_summ)
summary(fit_hh_size)

## Visualization

df_out_inf_all %>% 
  filter(r_beta == 0.25 & r_tau == 0.50 & r_omega == 0 & Inftot >=0)

df_out_inf_all %>% 
  filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.020 & time <= 70 & 
           n_hhsize == 5,
           n_exp_states == 3 & n_inf_states == 3,
         vax_prop == 0.9) #%>% View()

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


#### SMDM Figures ####
## Magnitude of epidemic peak
gg_epidemic_peak <- ggplot(df_out_inf_all_summ %>% 
         filter(r_beta == 0.25 & r_tau == 0.50 & r_omega == 0.000 & time <= 70 & 
                  n_exp_states == 3 & n_inf_states == 3 & n_hhsize %in% c(1, 3, 5) & 
                  eff_vax %in% c(0.9), vax_prop %in% c(0, 0.3, 0.6), level_npi != 0.8), 
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
gg_epidemic_peak
ggsave(plot = gg_epidemic_peak, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_peak_E3_I3.pdf", 
       width = 11, height = 8)
ggsave(plot = gg_epidemic_peak, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_peak_E3_I3.png", 
       width = 11, height = 8)
ggsave(plot = gg_epidemic_peak, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_peak_E3_I3.jpeg", 
       width = 11, height = 8)

## Epidemic curves
gg_epidemic_curve <- ggplot(df_out_inf_all %>% 
                             filter(r_beta == 0.25 & r_tau == 0.50 & r_omega == 0.000 & time <= 120 & 
                                      n_exp_states == 3 & n_inf_states == 3 & n_hhsize %in% c(1, 3, 5) & 
                                      eff_vax %in% c(0.5, 0.9), vax_prop %in% c(0, 0.3, 0.6), level_npi != 0.8), 
                           aes(x = time, y = Inftot/10e6, color = `Household size`, 
                               linetype = `Vaccine effectiveness`)) + # 
  geom_line(size = 1.1) +
  facet_grid(NPIeff ~ PropVax) +
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1.0)) +
  # scale_fill_viridis_d(option = "D", direction = -1) +
  scale_color_grey() +
  # scale_color_jcolors(palette = "rainbow") +
  # scale_fill_brewer(palette = "Spectral") +
  xlab("Time") +
  ylab("Infected population (% of total population)") +
  theme_bw(base_size = 16) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        legend.position = c("bottom"),
        # legend.position = c(.88,.8),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.key = element_blank())
gg_epidemic_curve
ggsave(plot = gg_epidemic_curve, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_curve_E3_I3.pdf", 
       width = 11, height = 8)
ggsave(plot = gg_epidemic_curve, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_curve_E3_I3.png", 
       width = 11, height = 8)
ggsave(plot = gg_epidemic_curve, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_curve_E3_I3.jpeg", 
       width = 11, height = 8)

gg_epidemic_curve_peak <- gg_epidemic_curve/gg_epidemic_peak
gg_epidemic_curve_peak <- gg_epidemic_curve_peak + plot_annotation(tag_levels = 'A')
gg_epidemic_curve_peak

ggsave(plot = gg_epidemic_curve_peak, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_curve_peak_E3_I3.pdf", 
       width = 11, height = 12)
ggsave(plot = gg_epidemic_curve_peak, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_curve_peak_E3_I3.png", 
       width = 11, height = 12)
ggsave(plot = gg_epidemic_curve_peak, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_curve_peak_E3_I3.jpeg", 
       width = 11, height = 12)




## Magnitude of epidemic peak
gg_epidemic_peak_E1_I1 <- ggplot(df_out_inf_all_summ %>% 
                             filter(r_beta == 0.25 & r_tau == 0.50 & r_omega == 0.000 & time <= 70 & 
                                      n_exp_states == 1 & n_inf_states == 1 & n_hhsize %in% c(1, 3, 5) & 
                                      eff_vax %in% c(0.9), vax_prop %in% c(0, 0.3, 0.6), level_npi != 0.8), 
                           aes(x = `Household size`, y = max_Inftot/10e6, fill = `Household size`)) + # linetype = as.factor(eff_vax))
  geom_col(color = NA) +
  facet_grid(NPIeff ~ PropVax) +
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1.0), limits = c(0, 0.06)) +
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
gg_epidemic_peak_E1_I1
ggsave(plot = gg_epidemic_peak_E1_I1, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_peak_E1_I1.pdf", 
       width = 11, height = 8)
ggsave(plot = gg_epidemic_peak_E1_I1, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_peak_E1_I1.png", 
       width = 11, height = 8)
ggsave(plot = gg_epidemic_peak_E1_I1, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_peak_E1_I1.jpeg", 
       width = 11, height = 8)

## Epidemic curves
gg_epidemic_curve_E1_I1 <- ggplot(df_out_inf_all %>% 
                              filter(r_beta == 0.25 & r_tau == 0.50 & r_omega == 0.000 & time <= 120 & 
                                       n_exp_states == 1 & n_inf_states == 1 & n_hhsize %in% c(1, 3, 5) & 
                                       eff_vax %in% c(0.5, 0.9), vax_prop %in% c(0, 0.3, 0.6), level_npi != 0.8), 
                            aes(x = time, y = Inftot/10e6, color = `Household size`, 
                                linetype = `Vaccine effectiveness`)) + # 
  geom_line(size = 1.1) +
  facet_grid(NPIeff ~ PropVax) +
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1.0), limits = c(0, 0.06)) +
  # scale_fill_viridis_d(option = "D", direction = -1) +
  scale_color_grey() +
  # scale_color_jcolors(palette = "rainbow") +
  # scale_fill_brewer(palette = "Spectral") +
  xlab("Time") +
  ylab("Infected population (% of total population)") +
  theme_bw(base_size = 16) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        # legend.position = c(""),
        legend.position = c(.88,.8),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.key = element_blank())
gg_epidemic_curve_E1_I1
ggsave(plot = gg_epidemic_curve_E1_I1, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_curve_E1_I1.pdf", 
       width = 11, height = 8)
ggsave(plot = gg_epidemic_curve_E1_I1, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_curve_E1_I1.png", 
       width = 11, height = 8)
ggsave(plot = gg_epidemic_curve_E1_I1, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_curve_E1_I1.jpeg", 
       width = 11, height = 8)

gg_epidemic_curve_peak_E1_I1 <- gg_epidemic_curve_E1_I1/gg_epidemic_peak_E1_I1
gg_epidemic_curve_peak_E1_I1 <- gg_epidemic_curve_peak_E1_I1 + plot_annotation(tag_levels = 'A')
gg_epidemic_curve_peak_E1_I1

ggsave(plot = gg_epidemic_curve_peak_E1_I1, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_curve_peak_E1_I1.pdf", 
       width = 11, height = 16)
ggsave(plot = gg_epidemic_curve_peak_E1_I1, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_curve_peak_E1_I1.png", 
       width = 11, height = 16)
ggsave(plot = gg_epidemic_curve_peak_E1_I1, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_curve_peak_E1_I1.jpeg", 
       width = 11, height = 16)

### Multiple E & I graphs

## Timing of epidemic peak
gg_epidemic_peak_time_E1_I1_E3_I3 <- ggplot(df_out_inf_all_summ %>% 
                                         filter(r_beta == 0.25 & r_tau == 0.50 & 
                                                  r_omega == 0.000 & time <= 70 & 
                                                  # n_exp_states %in% c(1, 3) & 
                                                  # n_inf_states %in% c(1, 3) & 
                                                  `Multicompartment structure` %in% c("E=1, I=1", "E=3, I=3") & 
                                                  n_hhsize %in% c(1, 5) & 
                                                  eff_vax %in% c(0.9) & 
                                                  vax_prop %in% c(0, 0.3, 0.6) & 
                                                  level_npi != 0.8), 
                                       aes(x = `Household size`, y = max_Inftot_time, 
                                           # fill = `Household size`, 
                                           fill = `Multicompartment structure`)) +
  geom_col(color = NA, position = "dodge") +
  facet_grid(NPIeff ~ PropVax) +
  scale_y_continuous() +
  # scale_fill_viridis_d(option = "D", direction = -1) +
  scale_fill_grey(start = 0.2, end = 0.6) +
  # scale_color_jcolors(palette = "rainbow") +
  # scale_fill_brewer(palette = "Spectral") +
  # xlab("Household size") +
  ylab("Timing of epidemic peak (days)") +
  guides(fill = guide_legend(nrow = 1)) +
  theme_bw(base_size = 16) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        # legend.position = c(""),
        legend.position = c(0.83, 0.8),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key = element_blank())
gg_epidemic_peak_time_E1_I1_E3_I3
ggsave(plot = gg_epidemic_peak_time_E1_I1_E3_I3, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_peak_time_E1_I1_E3_I3.pdf", 
       width = 11, height = 8)
ggsave(plot = gg_epidemic_peak_time_E1_I1_E3_I3, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_peak_time_E1_I1_E3_I3.png", 
       width = 11, height = 8)
ggsave(plot = gg_epidemic_peak_time_E1_I1_E3_I3, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_peak_time_E1_I1_E3_I3.jpeg", 
       width = 11, height = 8)

## Magnitude of epidemic peak
gg_epidemic_peak_E1_I1_E3_I3 <- ggplot(df_out_inf_all_summ %>% 
                                   filter(r_beta == 0.25 & r_tau == 0.50 & 
                                          r_omega == 0.000 & time <= 70 & 
                                          # n_exp_states %in% c(1, 3) & 
                                          # n_inf_states %in% c(1, 3) & 
                                          `Multicompartment structure` %in% c("E=1, I=1", "E=3, I=3") & 
                                          n_hhsize %in% c(1, 5) & 
                                          eff_vax %in% c(0.9) & 
                                          vax_prop %in% c(0, 0.3, 0.6) & 
                                          level_npi != 0.8), 
                                 aes(x = `Household size`, y = max_Inftot/10e6, 
                                     # fill = `Household size`, 
                                     fill = `Multicompartment structure`)) +
  geom_col(color = NA, position = "dodge") +
  facet_grid(NPIeff ~ PropVax) +
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1.0), limits = c(0, 0.06)) +
  # scale_fill_viridis_d(option = "D", direction = -1) +
  scale_fill_grey(start = 0.2, end = 0.6) +
  # scale_color_jcolors(palette = "rainbow") +
  # scale_fill_brewer(palette = "Spectral") +
  # xlab("Household size") +
  ylab("Magnitude of epidemic peak (% of total population)") +
  guides(fill = guide_legend(nrow = 1)) +
  theme_bw(base_size = 16) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        # legend.position = c(""),
        legend.position = c(0.83, 0.8),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key = element_blank())
gg_epidemic_peak_E1_I1_E3_I3
ggsave(plot = gg_epidemic_peak_E1_I1_E3_I3, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_peak_E1_I1_E3_I3.pdf", 
       width = 11, height = 8)
ggsave(plot = gg_epidemic_peak_E1_I1_E3_I3, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_peak_E1_I1_E3_I3.png", 
       width = 11, height = 8)
ggsave(plot = gg_epidemic_peak_E1_I1_E3_I3, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_peak_E1_I1_E3_I3.jpeg", 
       width = 11, height = 8)

## Epidemic curves
gg_epidemic_curve_E1_I1_E3_I3 <- ggplot(df_out_inf_all %>% 
                                    filter(r_beta == 0.25 & r_tau == 0.50 & 
                                           r_omega == 0.000 & time <= 120 & 
                                           # n_exp_states == 1 & 
                                           # n_inf_states == 1 & 
                                           `Multicompartment structure` %in% c("E=1, I=1", "E=3, I=3") & 
                                           n_hhsize %in% c(1, 5) & 
                                           eff_vax == 0.9 & # %in% c(0.5, 0.9) & 
                                           vax_prop %in% c(0, 0.3, 0.6) & 
                                           level_npi != 0.8), 
                                  aes(x = time, y = Inftot/10e6, color = `Household size`, 
                                      linetype = `Multicompartment structure`)) + # 
  geom_line(size = 1.1) +
  facet_grid(NPIeff ~ PropVax) +
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1.0), limits = c(0, 0.06)) +
  # scale_fill_viridis_d(option = "D", direction = -1) +
  scale_color_grey(start = 0.2, end = 0.6) +
  # scale_color_jcolors(palette = "rainbow") +
  # scale_fill_brewer(palette = "Spectral") +
  xlab("Time") +
  ylab("Infected population (% of total population)") +
  guides(color = guide_legend(nrow = 1), 
         linetype = guide_legend(nrow = 1)) +
  theme_bw(base_size = 16) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        # legend.position = c(""),
        legend.position = c(0.83, 0.84),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key = element_blank())
gg_epidemic_curve_E1_I1_E3_I3
ggsave(plot = gg_epidemic_curve_E1_I1_E3_I3, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_curve_E1_I1_E3_I3.pdf", 
       width = 11, height = 8)
ggsave(plot = gg_epidemic_curve_E1_I1_E3_I3, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_curve_E1_I1_E3_I3.png", 
       width = 11, height = 8)
ggsave(plot = gg_epidemic_curve_E1_I1_E3_I3, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_curve_E1_I1_E3_I3.jpeg", 
       width = 11, height = 8)


gg_epidemic_curve_peak_E1_I1_E3_I3 <- gg_epidemic_curve_E1_I1_E3_I3/gg_epidemic_peak_E1_I1_E3_I3
gg_epidemic_curve_peak_E1_I1_E3_I3 <- gg_epidemic_curve_peak_E1_I1_E3_I3 + plot_annotation(tag_levels = 'A')
gg_epidemic_curve_peak_E1_I1_E3_I3
ggsave(plot = gg_epidemic_curve_peak_E1_I1_E3_I3, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_curve_peak_E1_I1_E3_I3.pdf", 
       width = 11, height = 16)
ggsave(plot = gg_epidemic_curve_peak_E1_I1_E3_I3, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_curve_peak_E1_I1_E3_I3.png", 
       width = 11, height = 16)
ggsave(plot = gg_epidemic_curve_peak_E1_I1_E3_I3, 
       filename = "figs/SMDM_household_communiy_MC_SEIR_NPI_VAX_epidemic_curve_peak_E1_I1_E3_I3.jpeg", 
       width = 11, height = 16)
