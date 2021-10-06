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

## Figure 1: Model schematic ----

## Epidemic curves ----

### Figure 2: Epidemic Curves Natural History ----

max_Inftot_time_E1_I1 <- df_out_inf_all_summ %>% 
  filter(r_beta == 0.25 & r_tau == 0.50 & 
           r_omega == 0.000 & time <= 70 & 
           `Multicompartment structure` %in% c("E=1, I=1") & 
           n_hhsize %in% 1 & 
           eff_vax == 1 & 
           vax_prop == 0 & 
           level_npi == 1) %>%
  select(max_Inftot_time)

gg_epidemic_curve_nathist_E1_I1_E3_I3 <- ggplot(df_out_inf_all %>% 
                                                  filter(r_beta == 0.25 & r_tau == 0.50 & 
                                                           r_omega == 0.000 & time <= 70 & 
                                                           `Multicompartment structure` %in% c("E=1, I=1", "E=3, I=3", "E=3, I=1", "E=1, I=3") & 
                                                           n_hhsize %in% c(1, 3) & 
                                                           eff_vax == 1 & 
                                                           vax_prop == 0 & 
                                                           level_npi == 1), 
                                                aes(x = time, y = Inftot/10e6, color = `Household size`, 
                                                    linetype = `Multicompartment structure`)) + # 
  geom_line(size = 1.1) +
  geom_vline(xintercept = as.numeric(max_Inftot_time_E1_I1)) +
  facet_wrap(Esize ~ Isize) +
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
