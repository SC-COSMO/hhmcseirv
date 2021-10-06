#### Analyze epidemic outputs ####

rm(list = ls())

#### General setup ----
### Load packages ----
library(dplyr)
library(ggplot2)
library(patchwork)
library(stargazer)
# install.packages("remotes")
remotes::install_github("coolbutuseless/ggpattern")
source("R/02_decision_model_functions.R")

### Define color palette ----
jet.colors          <- colorRampPalette(c("black", "#00007F", 
                                          "blue", "#007FFF",
                                          "cyan", "#7FFF7F", 
                                          "yellow", "#FF7F00",
                                          "red", "#7F0000"))
jet.colors_mod      <- colorRampPalette(c("darkblue", "#03002e", 
                                          "blue", "#007FFF", 
                                          "cyan", "#7FFF7F", 
                                          "yellow", "#FF7F00", 
                                          "red", "#7F0000"))
color_map_mod       <- jet.colors_mod(120)
color_map           <- jet.colors(120)


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

## Calculate  differential outcomes compared to reference category (Group by parameter set)
# Create a data.frame with only reference category (ie, hhsize = 1)
df_out_inf_all_summ_ref <- df_out_inf_all_summ %>%
  filter(ref_category == 1) %>%
  rename(max_Inftot_ref = max_Inftot,
         max_Inftot_time_ref = max_Inftot_time,
         CumInfTot_ref = CumInfTot) %>%
  select(r_beta, r_tau, r_omega, 
         max_Inftot_ref, max_Inftot_time_ref, CumInfTot_ref)
# Compute differential outcomes compared to reference category (Group by parameter set)
df_out_inf_all_summ <- df_out_inf_all_summ %>%
  left_join(df_out_inf_all_summ_ref) %>%
  mutate(max_Inftot_diff = max_Inftot - max_Inftot_ref,
         max_Inftot_time_diff = max_Inftot_time - max_Inftot_time_ref,
         CumInfTot_diff = CumInfTot - CumInfTot_ref,
         # Percentage changes
         max_Inftot_diff_perc = (max_Inftot_diff/max_Inftot_ref)*100,
         max_Inftot_time_diff_perc = (max_Inftot_time_diff/max_Inftot_time_ref)*100,
         CumInfTot_diff_perc = (CumInfTot_diff/CumInfTot_ref)*100)

# hist(df_out_inf_all_summ$max_Inftot)
# hist(df_out_inf_all_summ$max_Inftot_time)
# hist(df_out_inf_all_summ$IDX100_time)
# hist(df_out_inf_all_summ$IDX500_time)
# hist(df_out_inf_all_summ$p05_Inftot_time, breaks = 15)
# hist(df_out_inf_all_summ$p10_Inftot_time, breaks = 15)
# hist(df_out_inf_all_summ$p25_Inftot_time, breaks = 15)
# hist(df_out_inf_all_summ$p50_Inftot_time, breaks = 15)

### Meta regression on category-specific outcomes
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

stargazer(fit_hh_peak_date, fit_hh_peak, fit_hh_size, 
          type = "text", out = "tables/03-metaregression_nathist_nointeractions.txt")
stargazer(fit_hh_peak_date, fit_hh_peak, fit_hh_size, 
          type = "latex", out = "tables/03-metaregression_nathist_nointeractions.tex")

## With interactions (for appendix)
fit_hh_int_peak_date <- lm(log(max_Inftot_time) ~ n_exp_states*n_inf_states + 
                         n_exp_states*n_hhsize + n_inf_states*n_hhsize +
                         r_tau + r_beta + r_omega, 
                       data = df_out_inf_all_summ)
summary(fit_hh_int_peak_date)

fit_hh_int_peak <- lm(log(max_Inftot) ~ n_exp_states*n_inf_states + 
                    n_exp_states*n_hhsize + n_inf_states*n_hhsize + 
                    n_hhsize + r_tau + r_beta + r_omega,
                  data = df_out_inf_all_summ)
summary(fit_hh_int_peak)

fit_hh_int_size <- lm(log(CumInfTot) ~ n_exp_states*n_inf_states + 
                    n_exp_states*n_hhsize + n_inf_states*n_hhsize +
                    n_hhsize + r_tau + r_beta + r_omega, 
                  data = df_out_inf_all_summ)
summary(fit_hh_int_size)

stargazer(fit_hh_int_peak_date, fit_hh_int_peak, fit_hh_int_size, 
          type = "text", 
          out = "tables/03-metaregression_nathist_interactions_appendix.txt")
stargazer(fit_hh_int_peak_date, fit_hh_int_peak, fit_hh_int_size, 
          type = "latex", 
          out = "tables/03-metaregression_nathist_interactions_appendix.tex")


### Meta regression on outcomes comparing to reference category 
## Without interactions
fit_hh_peak_date_diff <- lm(max_Inftot_time_diff ~ n_exp_states + 
                                   n_inf_states + 
                                   n_hhsize + 
                                   r_tau + r_beta + r_omega, 
                                 data = df_out_inf_all_summ)
summary(fit_hh_peak_date_diff)

fit_hh_peak_diff <- lm(max_Inftot_diff ~ n_exp_states + n_inf_states + 
                              n_hhsize*r_beta + 
                              r_tau + r_beta + r_omega,
                            data = df_out_inf_all_summ)
summary(fit_hh_peak_diff)

fit_hh_size_diff_perc <- lm(CumInfTot_diff_perc ~ n_exp_states + n_inf_states +
                              n_hhsize + r_tau + r_beta + r_omega, 
                            data = df_out_inf_all_summ)
summary(fit_hh_size_diff_perc)

stargazer(fit_hh_peak_date, fit_hh_peak, fit_hh_size, 
          type = "text", out = "tables/03-metaregression_nathist_nointeractions.txt")
stargazer(fit_hh_peak_date, fit_hh_peak, fit_hh_size, 
          type = "latex", out = "tables/03-metaregression_nathist_nointeractions.tex")

### Visualization
## Epidemic curves
gg_epidemic_curve_nathist <- ggplot(df_out_inf_all %>% 
                                          filter(r_beta == 0.25 & r_tau == 0.50 & 
                                                   r_omega == 0.000 & time <= 60 & 
                                                   # n_exp_states == 1 & 
                                                   # n_inf_states == 1 & 
                                                   # `Multicompartment structure` %in% c("E=1, I=1", "E=3, I=3") & 
                                                   n_hhsize < 7 &
                                                   level_npi == 1.0), 
                                        aes(x = time, y = Inftot/10e6, 
                                            color = `Household size`)) + # 
  geom_line(size = 1.05) +
  facet_grid(Esize ~ Isize) +
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1.0), limits = c(0, 0.06)) +
  scale_color_viridis_d(option = "C", direction = -1) +
  # scale_color_grey(start = 0.2, end = 0.6) +
  # scale_color_gradientn(colours = rev(color_map_mod)) +
  # scale_color_jcolors(palette = "rainbow") +
  # scale_fill_brewer(palette = "Spectral") +
  xlab("Time") +
  ylab("Infected population (% of total population)") +
  guides(color = guide_legend(title = "Household size", ncol = 1,
                              reverse = TRUE), 
         linetype = guide_legend(ncol = 1)) +
  theme_bw(base_size = 16) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        # legend.position = c(""),
        legend.position = c(0.94, 0.18),
        # legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.key = element_blank())
gg_epidemic_curve_nathist
ggsave(plot = gg_epidemic_curve_nathist, 
       filename = "figs/03_household_communiy_MC_SEIR_NATHIST_epidemic_curve.pdf", 
       width = 12, height = 8)
ggsave(plot = gg_epidemic_curve_nathist, 
       filename = "figs/03_household_communiy_MC_SEIR_NATHIST_epidemic_curve.pdf.png", 
       width = 12, height = 8)
ggsave(plot = gg_epidemic_curve_nathist, 
       filename = "figs/03_household_communiy_MC_SEIR_NATHIST_epidemic_curve.pdf.jpeg", 
       width = 12, height = 8)

## Epidemic measures vs reference category (ie, hhsize = 1, #e = 1, #i = 1)
ggplot(df_out_inf_all_summ %>%
         filter(r_beta == 0.25 & r_tau == 0.50 & 
                  r_omega == 0.000 & time <= 60 & 
                  n_hhsize < 7 &
                  level_npi == 1.0), 
       aes(x = `Household size`, 
           y = max_Inftot_diff_perc, 
           color = `Household size`)) + 
  geom_point() + 
  ylab("Percent difference in maximum number of infections") + 
  facet_grid(Esize ~ Isize) +
  theme(legend.position = "")
ggplot(df_out_inf_all_summ %>%
         filter(r_beta == 0.25 & r_tau == 0.50 & 
                  r_omega == 0.000 & time <= 60 & 
                  n_hhsize < 7 &
                  level_npi == 1.0), 
       aes(x = `Household size`, 
           y = max_Inftot_time_diff, 
           color = `Household size`)) + 
  geom_point() + 
  ylab("Difference in timing of maximum number of infections") + 
  facet_grid(Esize ~ Isize) +
  theme(legend.position = "")
ggplot(df_out_inf_all_summ %>%
         filter(r_beta == 0.25 & r_tau == 0.50 & 
                  r_omega == 0.000 & time <= 60 & 
                  n_hhsize < 7 &
                  level_npi == 1.0), 
       aes(x = `Household size`, 
           y = CumInfTot_diff, 
           color = `Household size`)) + 
  geom_point() + 
  ylab("Difference in total infections") + 
  facet_grid(Esize ~ Isize) +
  theme(legend.position = "")


# Outcomes by r_beta values.

ggplot(df_out_inf_all_summ %>%
         filter(r_tau == 0.40 & # r_beta == 0.25 & 
                r_omega == 0.000 & time <= 60 & 
                n_hhsize < 7 &
                level_npi == 1.0), 
       aes(x = `Household size`, 
           y = max_Inftot_diff_perc, 
           color = as.factor(r_beta))) + 
  geom_point() + 
  ylab("Percent difference in maximum number of infections") + 
  guides(color = guide_legend(title = "r_beta")) + 
  facet_grid(Esize ~ Isize) 
  # theme(legend.position = "")

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

## Calculate  differential outcomes compared to reference category (Group by parameter set)
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

## Calculate bias for all models compared bias to reference category (Group by parameter set)
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


## Compute differential outcomes compared to reference category
# df_out_inf_all_summ$max_Inftot_ref <- df_out_inf_all_summ$max_Inftot - 
#   df_out_inf_all_summ$max_Inftot[df_out_inf_all_summ$ref_category==1]
# df_out_inf_all_summ$max_Inftot_time_ref <- df_out_inf_all_summ$max_Inftot_time - 
#   df_out_inf_all_summ$max_Inftot_time[df_out_inf_all_summ$ref_category==1]
# df_out_inf_all_summ$CumInfTot_ref <- df_out_inf_all_summ$CumInfTot - 
#   df_out_inf_all_summ$CumInfTot[df_out_inf_all_summ$ref_category==1]


# Visualization ----

df_out_inf_all %>% 
  filter(r_beta == 0.25 & r_tau == 0.50 & r_omega == 0 & Inftot >=0)

df_out_inf_all %>% 
  filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.020 & time <= 70 & 
           n_hhsize == 5,
         n_exp_states == 3 & n_inf_states == 3,
         vax_prop == 0.9) #%>% View()

ggplot(df_out_inf_all %>% 
         filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.020 & time <= 70 & level_npi == 1), 
       aes(x = time, y = Inftot, color = `Household size`)) +
  geom_line(size = 1.3) +
  facet_grid(Esize ~ Isize) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom",
        strip.background = element_rect(colour="white", fill="lightgray"),
        # legend.position=c(.88,.3),
        legend.key = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin=margin(-10,-10,-10,-10))


# ggplot(df_out_inf_all %>% 
#          filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.000 & time <= 70 & 
#                   n_exp_states == 3 & n_inf_states == 3 & n_hhsize %in% c(1, 3, 5) & 
#                   eff_vax %in% c(0.9)), 
#        aes(x = time, y = Inftot, color = `Household size`)) + # linetype = as.factor(eff_vax))
#   geom_line(size = 1.1) +
#   facet_grid(NPIeff ~ PropVax) +
#   # scale_color_viridis_d(option = "C", direction = -1) +
#   theme_bw(base_size = 16) +
#   theme(legend.position = "bottom",
#         strip.background = element_rect(colour="white", fill="lightgray"),
#         # legend.position=c(.88,.3),
#         legend.key = element_blank(),
#         legend.margin = margin(0, 0, 0, 0),
#         legend.box.margin=margin(-10,-10,-10,-10))

#### Differential effect on control measures' effects
## TO-DO: https://stackoverflow.com/questions/28853786/how-do-i-plot-charts-with-nested-categories-axes/28868462#28868462
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

# ggplot(df_out_inf_all %>%
#          filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.020 & time <= 70 &
#                   n_exp_states == 3 & n_inf_states == 3 &
#                   eff_vax == 0.5),
#        aes(x = time, y = as.factor(eff_vax), fill = Inftot)) +
#   geom_tile() +
#   facet_grid((1-level_npi) + n_hhsize ~ vax_prop) +
#   theme_bw(base_size = 16) +
#   theme(legend.position = "bottom")

# List of to-dos:
# - add computation time
# - Add Incident infections equation and capture that on the output
# HH vs non-HH transmission.


#### Differential effect on control measures' effects BIAS
### Maximum infections
ggplot(df_out_inf_all_control_summ %>% 
         filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.000 & time <= 70 & 
                  n_hhsize %in% c(1, 3, 5) &
                  eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi != 1), 
       aes(x = `Household size`, y = max_Inftot_diff_bias_rel, 
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

### Time of maximum infections
ggplot(df_out_inf_all_control_summ %>% 
         filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.000 & time <= 70 & 
                  n_hhsize %in% c(1, 3, 5) &
                  eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi != 1), 
       aes(x = `Household size`, y = max_Inftot_time_diff_bias_rel, 
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

### Cumulative number of prevalent infections
ggplot(df_out_inf_all_control_summ %>% 
         filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.000 & time <= 70 & 
                  n_hhsize %in% c(1, 3, 5) &
                  eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi != 1), 
       aes(x = `Household size`, y = CumInfTot_diff_bias_rel, 
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
