#### Analyze epidemic outputs ####

rm(list = ls())

# General setup ----
## Load packages ----
library(dplyr)
library(ggplot2)
library(patchwork)
library(stargazer)
# install.packages("remotes")
# remotes::install_github("coolbutuseless/ggpattern")
# library(ggpattern) # https://coolbutuseless.github.io/2020/04/01/introducing-ggpattern-pattern-fills-for-ggplot/
source("R/02_decision_model_functions.R")

# Natural history outputs ----
## Load data ----
load(file = "output/df_output_doe_mc_seirv_all_nathist.RData")

## Rename variables for pretty plotting format ----
df_out_inf_all$Esize <- paste0("# of E compartments = ", df_out_inf_all$n_exp_states)
df_out_inf_all$Isize <- paste0("# of I compartments = ", df_out_inf_all$n_inf_states)
df_out_inf_all$`Household size` <- ordered(df_out_inf_all$n_hhsize, unique(df_out_inf_all$n_hhsize))
## df_out_inf_all$n_hhsize <- ordered(df_out_inf_all$n_hhsize)
df_out_inf_all$`Vaccine effectiveness` <- scales::percent(df_out_inf_all$eff_vax)
df_out_inf_all$PropVax <- paste0("Proportion vaccinated = ", scales::percent(df_out_inf_all$vax_prop))
df_out_inf_all$NPIeff <- paste0("NPI effectiveness = ", scales::percent(1-df_out_inf_all$level_npi))
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

## Epidemic curves ----

### Full Figure ----
gg_epidemic_curve_nathist_E1_I1_E3_I3 <- ggplot(df_out_inf_all %>% 
                                          filter(r_beta == 0.25 & r_tau == 0.50 & 
                                                   r_omega == 0.000 & time <= 70 & 
                                                   `Multicompartment structure` %in% c("E=1, I=1", "E=3, I=3") & 
                                                   n_hhsize %in% c(1, 3) & 
                                                   eff_vax == 1 & 
                                                   vax_prop == 0 & 
                                                   level_npi == 1), 
                                        aes(x = time, y = Inftot/10e6, color = `Household size`, 
                                            linetype = `Multicompartment structure`)) + # 
  geom_line(size = 1.1) +
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1.0), limits = c(0, 0.06)) +
  scale_color_grey(start = 0.2, end = 0.6) +
  xlab("Time") +
  ylab("Infected population (% of total population)") +
  guides(color = guide_legend(nrow = 1), 
         linetype = guide_legend(nrow = 1)) +
  theme_bw(base_size = 20) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        # legend.position = c(""),
        legend.position = c(0.77, 0.75),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key = element_blank())
gg_epidemic_curve_nathist_E1_I1_E3_I3
ggsave(plot = gg_epidemic_curve_nathist_E1_I1_E3_I3, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_epidemic_curve_nathist_E1_I1_E3_I3_all.pdf", 
       width = 8, height = 6)
ggsave(plot = gg_epidemic_curve_nathist_E1_I1_E3_I3, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_epidemic_curve_nathist_E1_I1_E3_I3_all.png", 
       width = 11, height = 8)
ggsave(plot = gg_epidemic_curve_nathist_E1_I1_E3_I3, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_epidemic_curve_nathist_E1_I1_E3_I3_all.jpeg", 
       width = 11, height = 8)

### HH = 1, E = 1, I = 1 ----
gg_epidemic_curve_nathist_E1_I1_hh1 <- ggplot(df_out_inf_all %>% 
                                                  filter(r_beta == 0.25 & r_tau == 0.50 & 
                                                           r_omega == 0.000 & time <= 70 & 
                                                           `Multicompartment structure` %in% c("E=1, I=1") & 
                                                           n_hhsize %in% c(1) & 
                                                           eff_vax == 1 & 
                                                           vax_prop == 0 & 
                                                           level_npi == 1), 
                                                aes(x = time, y = Inftot/10e6, color = `Household size`, 
                                                    linetype = `Multicompartment structure`)) + # 
  geom_line(size = 1.1) +
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1.0), limits = c(0, 0.06)) +
  scale_color_grey(start = 0.2, end = 0.6) +
  xlab("Time") +
  ylab("Infected population (% of total population)") +
  guides(color = guide_legend(nrow = 1, order = 1), 
         linetype = guide_legend(nrow = 1, order = 2)) +
  theme_bw(base_size = 20) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        # legend.position = c(""),
        legend.position = c(0.77, 0.75),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key = element_blank())
gg_epidemic_curve_nathist_E1_I1_hh1
ggsave(plot = gg_epidemic_curve_nathist_E1_I1_hh1, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_epidemic_curve_nathist_E1_I1_E3_I3_1ofX.pdf", 
       width = 8, height = 6)
ggsave(plot = gg_epidemic_curve_nathist_E1_I1_hh1, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_epidemic_curve_nathist_E1_I1_E3_I3_1ofX.png", 
       width = 11, height = 8)
ggsave(plot = gg_epidemic_curve_nathist_E1_I1_hh1, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_epidemic_curve_nathist_E1_I1_E3_I3_1ofX.jpeg", 
       width = 11, height = 8)

### HH = 1,3, E = 1, I = 1 ----
gg_epidemic_curve_nathist_E1_I1_hh1_3 <- ggplot(df_out_inf_all %>% 
                                                filter(r_beta == 0.25 & r_tau == 0.50 & 
                                                         r_omega == 0.000 & time <= 70 & 
                                                         `Multicompartment structure` %in% c("E=1, I=1") & 
                                                         n_hhsize %in% c(1, 3) & 
                                                         eff_vax == 1 & 
                                                         vax_prop == 0 & 
                                                         level_npi == 1), 
                                              aes(x = time, y = Inftot/10e6, color = `Household size`, 
                                                  linetype = `Multicompartment structure`)) + # 
  geom_line(size = 1.1) +
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1.0), limits = c(0, 0.06)) +
  scale_color_grey(start = 0.2, end = 0.6) +
  xlab("Time") +
  ylab("Infected population (% of total population)") +
  guides(color = guide_legend(nrow = 1, order = 1), 
         linetype = guide_legend(nrow = 1, order = 2)) +
  theme_bw(base_size = 20) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        # legend.position = c(""),
        legend.position = c(0.77, 0.75),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key = element_blank())
gg_epidemic_curve_nathist_E1_I1_hh1_3
ggsave(plot = gg_epidemic_curve_nathist_E1_I1_hh1_3, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_epidemic_curve_nathist_E1_I1_E3_I3_2ofX.pdf", 
       width = 8, height = 6)
ggsave(plot = gg_epidemic_curve_nathist_E1_I1_hh1_3, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_epidemic_curve_nathist_E1_I1_E3_I3_2ofX.png", 
       width = 11, height = 8)
ggsave(plot = gg_epidemic_curve_nathist_E1_I1_hh1_3, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_epidemic_curve_nathist_E1_I1_E3_I3_2ofX.jpeg", 
       width = 11, height = 8) 

### HH = 1, E = 1,3, I = 1,3 ----
gg_epidemic_curve_nathist_E1_I1_E3_I3 <- ggplot(df_out_inf_all %>% 
                                                filter(r_beta == 0.25 & r_tau == 0.50 & 
                                                         r_omega == 0.000 & time <= 70 & 
                                                         `Multicompartment structure` %in% c("E=1, I=1", "E=3, I=3") & 
                                                         n_hhsize %in% c(1) & 
                                                         eff_vax == 1 & 
                                                         vax_prop == 0 & 
                                                         level_npi == 1), 
                                              aes(x = time, y = Inftot/10e6, color = `Household size`, 
                                                  linetype = `Multicompartment structure`)) + # 
  geom_line(size = 1.1) +
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1.0), limits = c(0, 0.06)) +
  scale_color_grey(start = 0.2, end = 0.6) +
  xlab("Time") +
  ylab("Infected population (% of total population)") +
  guides(color = guide_legend(nrow = 1, order = 1), 
         linetype = guide_legend(nrow = 1, order = 2)) +
  theme_bw(base_size = 20) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        # legend.position = c(""),
        legend.position = c(0.77, 0.75),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key = element_blank())
gg_epidemic_curve_nathist_E1_I1_E3_I3
ggsave(plot = gg_epidemic_curve_nathist_E1_I1_E3_I3, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_epidemic_curve_nathist_E1_I1_E3_I3_3ofX.pdf", 
       width = 8, height = 6)
ggsave(plot = gg_epidemic_curve_nathist_E1_I1_E3_I3, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_epidemic_curve_nathist_E1_I1_E3_I3_3ofX.png", 
       width = 8, height = 6)
ggsave(plot = gg_epidemic_curve_nathist_E1_I1_E3_I3, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_epidemic_curve_nathist_E1_I1_E3_I3_3ofX.jpeg", 
       width = 8, height = 6)

## Meta regression on category-specific outcomes ----
## Without interactions
fit_hh_peak_date <- lm(log(max_Inftot_time) ~ n_exp_states + n_inf_states + 
                         n_hhsize + 
                         r_tau + r_beta + r_omega, 
                       data = df_out_inf_all_summ)
summary(fit_hh_peak_date)

fit_hh_peak <- lm(log(max_Inftot) ~ n_exp_states + n_inf_states + 
                    n_hhsize + 
                    r_tau + r_beta + r_omega,
                  data = df_out_inf_all_summ)
summary(fit_hh_peak)

fit_hh_size <- lm(log(CumInfTot) ~ n_exp_states + n_inf_states + 
                    n_hhsize + r_tau + r_beta + r_omega, 
                  data = df_out_inf_all_summ)
summary(fit_hh_size)

stargazer(fit_hh_peak_date, fit_hh_peak, fit_hh_size, 
          type = "text")

# Control measure (vaccination and NPI) outputs ----

## Load data ----
load(file = "output/df_output_doe_mc_seirv_all_control.RData")

## Wrangle data ----
# Rename variables for pretty plotting format
df_out_inf_all$Esize <- paste0("# of E compartments = ", df_out_inf_all$n_exp_states)
df_out_inf_all$Isize <- paste0("# of I compartments = ", df_out_inf_all$n_inf_states)
df_out_inf_all$`Household size` <- ordered(df_out_inf_all$n_hhsize, unique(df_out_inf_all$n_hhsize))
# df_out_inf_all$n_hhsize <- ordered(df_out_inf_all$n_hhsize)
df_out_inf_all$`Vaccine effectiveness` <- scales::percent(df_out_inf_all$eff_vax)
df_out_inf_all$PropVax <- paste0("Proportion vaccinated = ", scales::percent(df_out_inf_all$vax_prop))
df_out_inf_all$NPIeff <- paste0("NPI effectiveness = ", scales::percent(1-df_out_inf_all$level_npi))
df_out_inf_all$NPIeff_simple <- paste0(scales::percent(1-df_out_inf_all$level_npi))
df_out_inf_all$`Multicompartment structure` <- paste0("E=", 
                                                      df_out_inf_all$n_exp_states, 
                                                      ", I=", 
                                                      df_out_inf_all$n_inf_states)

## Summarize output ----
df_out_inf_all_control_summ <- df_out_inf_all %>%
  group_by(pid) %>% # ,  n_exp_states, n_inf_states
  mutate(# Find the t at which I(t) is at its max
    ref_category = ifelse(n_hhsize == 1, 
                          1, 0), # Define HH=1 as reference category
    ref_category_bias = ifelse(n_hhsize == 1 & n_exp_states == 1 & n_inf_states == 1.0, 
                               1, 0), # n_hhsize = 1, E = 1, I =1 as reference category to estimate bias of control measures' effect
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
    IDX500_time = max(time[which((Inftot-InfNoDX) <= 500 & time < max_Inftot_time)]),
    IDX100_time = max(time[which((Inftot-InfNoDX) <= 100 & time < max_Inftot_time)]),
    CumInfTot = sum(Inftot)) %>%
  # arrange(-ref_category_nathist) %>%
  slice_head() %>%
  ungroup()

## Calculate  differential outcomes compared to reference category (Group by parameter set) ----
# Create a data.frame with only reference category (ie, hhsize = 1)
df_out_inf_all_control_summ_nathist <- df_out_inf_all_control_summ %>%
  filter(ref_category_nathist == 1) %>%
  rename(max_Inftot_nathist = max_Inftot,
         max_Inftot_time_nathist = max_Inftot_time,
         CumInfTot_nathist = CumInfTot) %>%
  select(n_hhsize, n_exp_states, n_inf_states,
         r_beta, r_tau, r_omega, 
         max_Inftot_nathist, max_Inftot_time_nathist, CumInfTot_nathist)

# Compute differential outcomes compared to reference category (Group by parameter set)
df_out_inf_all_control_summ <- df_out_inf_all_control_summ %>%
  left_join(df_out_inf_all_control_summ_nathist) %>%
  mutate(max_Inftot_diff = max_Inftot - max_Inftot_nathist,
         max_Inftot_time_diff = max_Inftot_time - max_Inftot_time_nathist,
         CumInfTot_diff = CumInfTot - CumInfTot_nathist,
         # Percentage changes
         max_Inftot_diff_perc = (max_Inftot_diff/max_Inftot_nathist)*100,
         max_Inftot_time_diff_perc = (max_Inftot_time_diff/max_Inftot_time_nathist)*100,
         CumInfTot_diff_perc = (CumInfTot_diff/CumInfTot_nathist)*100)

## Calculate bias for all models compared bias to reference category (Group by parameter set) ----
# Create a data.frame with only reference category (ie, hhsize = 1, e =1, i = 1)
df_out_inf_all_control_summ_ref_bias <- df_out_inf_all_control_summ %>%
  filter(ref_category_bias == 1) %>%
  rename(max_Inftot_diff_ref_bias = max_Inftot_diff,
         max_Inftot_time_diff_ref_bias = max_Inftot_time_diff,
         CumInfTot_diff_ref_bias = CumInfTot_diff) %>%
  select(r_beta, r_tau, r_omega, 
         max_Inftot_diff_ref_bias, 
         max_Inftot_time_diff_ref_bias, 
         CumInfTot_diff_ref_bias, vax_prop, 
         level_npi, eff_vax)

# Compute bias in differential outcomes compared to bia reference category (Group by parameter set)
df_out_inf_all_control_summ <- df_out_inf_all_control_summ %>%
  left_join(df_out_inf_all_control_summ_ref_bias) %>%
  mutate(# Absolute bias
    max_Inftot_diff_bias_abs = max_Inftot_diff_ref_bias - max_Inftot_diff,
    max_Inftot_time_diff_bias_abs = max_Inftot_time_diff_ref_bias - max_Inftot_time_diff,
    CumInfTot_diff_bias_abs = CumInfTot_diff_ref_bias - CumInfTot_diff,
    # Relative bias
    max_Inftot_diff_bias_rel = 100 * max_Inftot_diff_bias_abs/max_Inftot_diff,
    max_Inftot_time_diff_bias_rel = 100 * max_Inftot_time_diff_bias_abs/max_Inftot_time_diff,
    CumInfTot_diff_bias_rel = 100 * CumInfTot_diff_bias_abs/CumInfTot_diff)

## Control measures' effects ----
# Format data for plotting
df_control_measures_E1_I1_E3_I3 <- df_out_inf_all_control_summ %>% 
  filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.000 & time <= 70 & 
           ((`Multicompartment structure` %in% c("E=1, I=1") & n_hhsize %in% c(1)) | 
              (`Multicompartment structure` %in% c("E=3, I=3") & n_hhsize %in% c(3))) & 
           eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi %in% c(0.4, 1.0)) %>%
  mutate(`Natural history structure`  = paste0(`Multicompartment structure`, " & HH = ", n_hhsize))

df_control_measures_E1_I1_E3_I3 %>% 
  select(`Multicompartment structure`, n_hhsize, level_npi, max_Inftot)

### Full Figure ----
gg_control_measures_E1_I1_E3_I3 <- ggplot(df_control_measures_E1_I1_E3_I3,
                                          aes(x = `Natural history structure`, y = max_Inftot,
                                              group = NPIeff, fill = NPIeff_simple)) +
  geom_col(width = 0.5, position = position_dodge(0.5)) +
  scale_fill_grey("NPI Effectiveness") +
  scale_y_continuous("Peak size (thousands)", breaks = seq(0, 600000, by = 100000),
                     labels = function(x) round(x/1000, digits = 0)) +
  scale_x_discrete() + 
  # xlab("Household size") +
  theme_bw(base_size = 20) +
  # coord_flip() +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        axis.title.x=element_blank(),
        legend.position = c(0.5, 0.3),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.key = element_blank())
gg_control_measures_E1_I1_E3_I3
ggsave(plot = gg_control_measures_E1_I1_E3_I3, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_E1_I1_E3_I3_all.pdf", 
       width = 12, height = 6)
ggsave(plot = gg_control_measures_E1_I1_E3_I3, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_E1_I1_E3_I3_all.png", 
       width = 12, height = 6)
ggsave(plot = gg_control_measures_E1_I1_E3_I3, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_E1_I1_E3_I3_all.jpeg", 
       width = 12, height = 6)

### HH = 1,3, E = 1,3, I = 1,3, NPI Eff = 0% ----
gg_control_measures_E1_I1_E3_I3_NPIeff0 <- ggplot(df_control_measures_E1_I1_E3_I3 %>%
                                                    filter(NPIeff_simple == "0%"),
                                          aes(x = `Natural history structure`, y = max_Inftot,
                                              group = NPIeff, fill = NPIeff_simple)) +
  geom_col(width = 0.5, position = position_dodge(0.5)) +
  scale_fill_grey("NPI Effectiveness") +
  scale_y_continuous("Peak size (thousands)", breaks = seq(0, 600000, by = 100000),
                     labels = function(x) round(x/1000, digits = 0)) +
  scale_x_discrete() + 
  # xlab("Household size") +
  ylab("Epidemic size (millions)") +
  theme_bw(base_size = 20) +
  # coord_flip() +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        axis.title.x=element_blank(),
        legend.position = c(0.5, 0.3),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.key = element_blank())
gg_control_measures_E1_I1_E3_I3_NPIeff0
ggsave(plot = gg_control_measures_E1_I1_E3_I3_NPIeff0, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_E1_I1_E3_I3_1ofX.pdf", 
       width = 12, height = 6)
ggsave(plot = gg_control_measures_E1_I1_E3_I3_NPIeff0, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_E1_I1_E3_I3_1ofX.png", 
       width = 12, height = 6)
ggsave(plot = gg_control_measures_E1_I1_E3_I3_NPIeff0, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_E1_I1_E3_I3_1ofX.jpeg", 
       width = 12, height = 6)

## Differential effect on control measures' effects ----
ggplot(df_out_inf_all_control_summ %>% 
         filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.020 & time <= 70 & 
                  n_hhsize %in% c(1, 3, 5) &
                  eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi != 1), 
       aes(x = `Household size`, y = -max_Inftot_diff, 
           group = Isize, color = Isize, fill = Esize)) + # linetype = as.factor(eff_vax))
  geom_bar(stat = "identity", position = position_dodge(0.8)) +
  facet_grid(NPIeff ~ Esize) +
  # facet_wrap(NPIeff ~ Esize, scales = "free") +
  scale_fill_grey() +
  # scale_y_continuous(labels = function(x)round(x/10e6, digits = 2)) +
  # scale_fill_viridis_d(option = "C", direction = -1) +
  # scale_color_jcolors(palette = "rainbow") +
  # scale_fill_brewer(palette = "Spectral") +
  # xlab("Household size") +
  # ylab("Cumulative infections (millions)") +
  theme_bw(base_size = 16) +
  # coord_flip() +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        # legend.position = c(""),
        # legend.position = c(.88,.8),
        # legend.position = "bottom",
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin=margin(-10,-10,-10,-10)
        legend.key = element_blank())

## Differential effect on control measures' effects BIAS ----
alphas <- c("# of I compartments = 1" = 0.4, 
            "# of I compartments = 2" = 0.7, 
            "# of I compartments = 3" = 1.0)
### Maximum infections ----
gg_peak_size_bias_rel <- ggplot(df_out_inf_all_control_summ %>% 
         filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.000 & time <= 70 & 
                  n_hhsize %in% c(1, 3, 5) &
                  eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi != 1), 
       aes(x = `Household size`, y = max_Inftot_diff_bias_rel, 
           group = Isize, alpha = Isize, fill = Esize)) + # linetype = as.factor(eff_vax))
  geom_bar(stat = "identity", position = position_dodge(0.8)) +
  facet_grid(NPIeff ~ Esize) +
  scale_fill_discrete(l = 50) +
  scale_alpha_manual(values = alphas) +
  ylab("Relative bias on peak size (%)") +
  theme_bw(base_size = 16) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        legend.position = c(.79, .86),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank())
gg_peak_size_bias_rel
ggsave(plot = gg_peak_size_bias_rel, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_peak_size_bias_rel.pdf", 
       width = 12, height = 8)
ggsave(plot = gg_peak_size_bias_rel, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_peak_size_bias_rel.png", 
       width = 12, height = 8)
ggsave(plot = gg_peak_size_bias_rel, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_peak_size_bias_rel.jpeg", 
       width = 12, height = 8)

### Time of maximum infections
gg_peak_time_bias_rel <- ggplot(df_out_inf_all_control_summ %>% 
         filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.000 & time <= 70 & 
                  n_hhsize %in% c(1, 3, 5) &
                  eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi != 1), 
       aes(x = `Household size`, y = max_Inftot_time_diff_bias_rel, 
           group = Isize, alpha = Isize, fill = Esize)) + # linetype = as.factor(eff_vax))
  geom_bar(stat = "identity", position = position_dodge(0.8)) +
  facet_grid(NPIeff ~ Esize) +
  scale_fill_discrete(l = 50) +
  scale_alpha_manual(values = alphas) +
  # xlab("Household size") +
  scale_y_continuous("Relative bias on peak time (%)", 
                     breaks = c(0, 250, 500, 750, 1000, 1250), 
                     labels = function(x) scales::comma(x)) +
  theme_bw(base_size = 16) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        legend.position = c(.79, .86),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank())
gg_peak_time_bias_rel
ggsave(plot = gg_peak_time_bias_rel, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_peak_time_bias_rel.pdf", 
       width = 12, height = 8)
ggsave(plot = gg_peak_time_bias_rel, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_peak_time_bias_rel.png", 
       width = 12, height = 8)
ggsave(plot = gg_peak_time_bias_rel, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_peak_time_bias_rel.jpeg", 
       width = 12, height = 8)

### Cumulative number of prevalent infections
gg_epidemic_size_bias_rel <- ggplot(df_out_inf_all_control_summ %>% 
                                  filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.000 & time <= 70 & 
                                           n_hhsize %in% c(1, 3, 5) &
                                           eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi != 1), 
                                aes(x = `Household size`, y = CumInfTot_diff_bias_rel, 
                                    group = Isize, alpha = Isize, fill = Esize)) + # linetype = as.factor(eff_vax))
  geom_bar(stat = "identity", position = position_dodge(0.8)) +
  facet_grid(NPIeff ~ Esize) +
  scale_fill_discrete(l = 50) +
  scale_alpha_manual(values = alphas) +
  # xlab("Household size") +
  scale_y_continuous("Relative bias on epidemic size (%)", 
                     breaks = seq(0, 700, by = 100), 
                     labels = function(x) scales::comma(x)) +
  theme_bw(base_size = 16) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        legend.position = c(.79, .86),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank())
gg_epidemic_size_bias_rel
ggsave(plot = gg_epidemic_size_bias_rel, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_epidemic_size_bias_rel.pdf", 
       width = 12, height = 8)
ggsave(plot = gg_epidemic_size_bias_rel, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_epidemic_size_bias_rel.png", 
       width = 12, height = 8)
ggsave(plot = gg_epidemic_size_bias_rel, 
       filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_epidemic_size_bias_rel.jpeg", 
       width = 12, height = 8)

## Meta regression on category-specific outcomes ----
df_control_measures_bias <- df_out_inf_all_control_summ %>% 
  filter(r_omega == 0.000 & 
           !((`Multicompartment structure` %in% c("E=1, I=1") & n_hhsize %in% c(1))) & 
           eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi %in% c(0.4)) 

## Without interactions
fit_hh_peak_time_bias_rel <- lm(max_Inftot_time_diff_bias_rel ~ n_exp_states + n_inf_states + 
                         n_hhsize + 
                         r_tau + r_beta + r_omega, 
                       data = df_control_measures_bias %>% filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel)

fit_hh_peak_size_bias_rel <- lm(max_Inftot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                    n_hhsize + 
                    r_tau + r_beta + r_omega,
                  data = df_control_measures_bias)
summary(fit_hh_peak_size_bias_rel)

fit_hh_epidemic_size <- lm(CumInfTot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                    n_hhsize + r_tau + r_beta + r_omega, 
                  data = df_control_measures_bias)
summary(fit_hh_epidemic_size)

stargazer(fit_hh_peak_time_bias_rel, 
          fit_hh_peak_size_bias_rel, 
          fit_hh_epidemic_size, 
          type = "text")

## Meta regression on category-specific outcomes and absolute value of the relative bias----
df_control_measures_bias <- df_out_inf_all_control_summ %>% 
  filter(r_omega == 0.000 & 
           !((`Multicompartment structure` %in% c("E=1, I=1") & n_hhsize %in% c(1))) & 
           eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi %in% c(0.4)) 

## Without interactions
fit_hh_peak_time_bias_rel_abs <- lm(abs(max_Inftot_time_diff_bias_rel) ~ n_exp_states + n_inf_states + 
                                  n_hhsize + 
                                  r_tau + r_beta , 
                                data = df_control_measures_bias %>% 
                                  filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel_abs)

fit_hh_peak_size_bias_rel_abs <- lm(abs(max_Inftot_diff_bias_rel) ~ n_exp_states + n_inf_states + 
                                  n_hhsize + 
                                  r_tau + r_beta + r_omega,
                                data = df_control_measures_bias)
summary(fit_hh_peak_size_bias_rel_abs)

fit_hh_epidemic_size_abs <- lm(abs(CumInfTot_diff_bias_rel) ~ n_exp_states + n_inf_states + 
                             n_hhsize + r_tau + r_beta + r_omega, 
                           data = df_control_measures_bias)
summary(fit_hh_epidemic_size_abs)

stargazer(fit_hh_peak_time_bias_rel_abs, 
          fit_hh_peak_size_bias_rel_abs, 
          fit_hh_epidemic_size_abs, 
          type = "text")



df_control_measures_bias %>% 
  filter(max_Inftot_time_diff_bias_rel != Inf) %>%
  hist(max_Inftot_diff_bias_rel)

hist(df_control_measures_bias$max_Inftot_time_diff_bias_rel)
hist(df_control_measures_bias$max_Inftot_diff_bias_rel)
hist(df_control_measures_bias$CumInfTot_diff_bias_rel)
