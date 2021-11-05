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

max_Inftot_time_E1_I1_hh1 <- df_out_inf_all_summ %>% 
  filter(r_beta == 0.25 & r_tau == 0.50 & 
           r_omega == 0.000 & time <= 70 & 
           `Multicompartment structure` %in% c("E=1, I=1") & 
           n_hhsize %in% 1 & 
           eff_vax == 1 & 
           vax_prop == 0 & 
           level_npi == 1) %>%
  select(max_Inftot_time, max_Inftot)
max_Inftot_time_E1_I1_hh1$max_Inftot <- max_Inftot_time_E1_I1_hh1$max_Inftot/10e6

max_Inftot_time_E1_I1_hh3 <- df_out_inf_all_summ %>% 
  filter(r_beta == 0.25 & r_tau == 0.50 & 
           r_omega == 0.000 & time <= 70 & 
           `Multicompartment structure` %in% c("E=1, I=1") & 
           n_hhsize %in% 3 & 
           eff_vax == 1 & 
           vax_prop == 0 & 
           level_npi == 1) %>%
  select(max_Inftot_time, max_Inftot)
max_Inftot_time_E1_I1_hh3$max_Inftot <- max_Inftot_time_E1_I1_hh3$max_Inftot/10e6

max_Inftot_time_E1_I1 <- df_out_inf_all_summ %>% 
  filter(r_beta == 0.25 & r_tau == 0.50 & 
           r_omega == 0.000 & time <= 70 & 
           `Multicompartment structure` %in% c("E=1, I=1", 
                                               "E=3, I=3", 
                                               "E=3, I=1", 
                                               "E=1, I=3") & 
           n_hhsize %in% c(1, 3) & 
           eff_vax == 1 & 
           vax_prop == 0 & 
           level_npi == 1) %>%
  mutate(x0 = c(rep(max_Inftot_time_E1_I1_hh1$max_Inftot_time, 4),
                rep(max_Inftot_time_E1_I1_hh3$max_Inftot_time, 4))) %>% 
  mutate(y0 = c(rep(max_Inftot_time_E1_I1_hh1$max_Inftot, 4),
                rep(max_Inftot_time_E1_I1_hh3$max_Inftot, 4))) %>% 
  mutate(origin = (x0 == max_Inftot_time),
         x0 = ifelse(origin, NA, x0),
         y0 = ifelse(origin, NA, y0)) %>%
  select(`Household size`,
         `Multicompartment structure`, 
         max_Inftot_time, 
         max_Inftot,
         origin,
         x0, y0)

max_Inftot_time_E1_I1$max_Inftot <- max_Inftot_time_E1_I1$max_Inftot/10e6

df_fig2 <- left_join(df_out_inf_all %>% 
                       filter(r_beta == 0.25 & r_tau == 0.50 & 
                                r_omega == 0.000 & time <= 60 & 
                                `Multicompartment structure` %in% c("E=1, I=1", 
                                                                    "E=3, I=3", 
                                                                    "E=3, I=1", 
                                                                    "E=1, I=3") & 
                                n_hhsize %in% c(1, 3) & 
                                eff_vax == 1 & 
                                vax_prop == 0 & 
                                level_npi == 1),
                     max_Inftot_time_E1_I1)

gg_epidemic_curve_nathist_E1_I1_E3_I3 <- ggplot(df_fig2, 
                                                aes(x = time, y = Inftot/10e6, color = `Household size`)) + # 
  geom_line(size = 1.1) +
  geom_segment(aes(x = max_Inftot_time_E1_I1_hh1$max_Inftot_time, 
                   xend = max_Inftot_time_E1_I1_hh1$max_Inftot_time, 
                   y = 0, 
                   yend = max_Inftot_time_E1_I1_hh1$max_Inftot),
               linetype = "dashed", color = "blue") +
  geom_segment(aes(x = max_Inftot_time_E1_I1_hh3$max_Inftot_time, 
                   xend = max_Inftot_time_E1_I1_hh3$max_Inftot_time, 
                   y = 0, 
                   yend = max_Inftot_time_E1_I1_hh3$max_Inftot),
               linetype = "dashed", color = "red") +
  geom_point(aes(x = max_Inftot_time_E1_I1_hh1$max_Inftot_time, 
                 y = max_Inftot_time_E1_I1_hh1$max_Inftot),
             color = "blue", size = 2.5) +
  geom_point(aes(x = max_Inftot_time_E1_I1_hh3$max_Inftot_time, 
                 y = max_Inftot_time_E1_I1_hh3$max_Inftot),
             color = "red", size = 2.5) +
  geom_segment(aes(x = x0, xend = max_Inftot_time,
                   y = y0, yend = max_Inftot,
                   color = `Household size`),
               linejoin = "mitre",
               linetype = "dashed", arrow = arrow(length = unit(0.15, "inches"),
                                                  type = "closed"), 
               show.legend = FALSE) +
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
gg_epidemic_curve_nathist_E1_I1_E3_I3

ggsave(plot = gg_epidemic_curve_nathist_E1_I1_E3_I3, 
       filename = "figs/Paper/Fig2_natural_history_curves.pdf", 
       width = 12, height = 8)

### Table 2: Meta regression on category-specific outcomes on Natural History ----
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

## With interactions
fit_hh_peak_date_int <- lm(log(max_Inftot_time) ~ n_exp_states*n_inf_states + 
                             n_exp_states*n_hhsize + 
                             n_inf_states*n_hhsize + 
                         r_tau + r_beta + r_omega, 
                       data = df_out_inf_all_summ)
summary(fit_hh_peak_date_int)

fit_hh_peak_int <- lm(log(max_Inftot) ~ n_exp_states*n_inf_states + 
                        n_exp_states*n_hhsize + 
                        n_inf_states*n_hhsize +
                    r_tau + r_beta + r_omega,
                  data = df_out_inf_all_summ)
summary(fit_hh_peak_int)

fit_hh_size_int <- lm(log(CumInfTot) ~ n_exp_states*n_inf_states + 
                        n_exp_states*n_hhsize + 
                        n_inf_states*n_hhsize + 
                        r_tau + r_beta + r_omega, 
                  data = df_out_inf_all_summ)
summary(fit_hh_size_int)

stargazer(fit_hh_peak_date, fit_hh_peak_date_int,
          fit_hh_peak, fit_hh_peak_int, 
          fit_hh_size, fit_hh_size_int, 
          type = "text", # latex 
          title = "Metaregression estimates",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          covariate.labels = c("E", "I", "Household size (HH)", "$\\tau$", "$\\beta$", "$\\omega$", "E*I", "E*HH", "I*HH"), 
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE) 
stargazer(fit_hh_peak_date, fit_hh_peak_date_int,
          fit_hh_peak, fit_hh_peak_int, 
          fit_hh_size, fit_hh_size_int,
          type = "latex",  
          out =  "output/Table2_metaregression_epidemic_curves.tex",
          title = "Metaregression estimates",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          covariate.labels = c("E", "I", "Household size (HH)", "$\\tau$", "$\\beta$", "$\\omega$", "E*I", "E*HH", "I*HH"), 
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE) 
# Control measure outputs ----
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
df_out_inf_all$NPIeff_labels <- ordered(df_out_inf_all$NPIeff,
                                        unique(df_out_inf_all$NPIeff), c("No NPI", "NPI20", "NPI"))
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


### Figure 3: Differential effect on control measures' effects BIAS ----
png(file = "figs/Paper/Fig3_exemplary_bias.png", width = 960, height = 960*0.8)
gg_control_measures_E1_I1_E3_I3 <- ggplot(df_control_measures_E1_I1_E3_I3,
                                          aes(x = `Natural history structure`, y = max_Inftot,
                                              group = NPIeff, fill = NPIeff_labels)) +
  geom_col(width = 0.5, position = position_dodge(0.5)) +
  # scale_fill_grey("Intervention") +
  scale_fill_manual("Intervention", values = c("No NPI" = "darkgreen", "NPI" = "darkred")) + 
  scale_y_continuous("Peak size (thousands)", breaks = seq(0, 600000, by = 100000),
                     labels = function(x) round(x/1000, digits = 0)) +
  scale_x_discrete() + 
  geom_bracket(y.position = 410*1000, label = "O[NH]", type = "expression", 
               xmin = 0.7, 
               xmax = 1.3, inherit.aes = FALSE, ) +
  geom_bracket(y.position = 580*1000, label = "O[HH]", type = "expression", 
               xmin = 1.7, 
               xmax = 2.3, inherit.aes = FALSE) +
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
## First pair of brackets
grid.brackets(x1 = 250, x2 = 250, 
              y1 = 385, y2 = 275, 
              lwd = 1)
grid.text(x = unit(282, "native"), y = unit(325, "native"),
          label = expression(paste(Delta, O[NH])), hjust = 0, vjust=0)
## Second pair of brackets
grid.brackets(x1 = 540, x2 = 540, 
              y1 = 525, y2 = 445, 
              lwd = 1)
grid.text(x = unit(572, "native"), y = unit(480, "native"),
          label = expression(paste(Delta, O[HH])), hjust = 0, vjust=0)

grid.text(x = unit(100, "native"), y = unit(500, "native"),
          label =  expression(paste("rBias = ", 
                                    frac(paste(Delta, O[NH]) - paste(Delta, O[HH]), 
                                         paste(Delta, O[HH])
                                    )%*%100, " = 30%")
          ), 
          hjust = 0, vjust=0)
dev.off()

#### Maximum infections ----
alphas <- c("# of I compartments = 1" = 0.4, 
            "# of I compartments = 2" = 0.7, 
            "# of I compartments = 3" = 1.0)
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
  theme_bw(base_size = 20) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        legend.position = c(.84, .86),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank())
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
          font.label=list(color = "black", size = 16))

ggsave(plot = fig3, 
       filename = "figs/Paper/Fig3.pdf", 
       width = 24, height = 10)

### Appendix figures: Differential effect on control measures' effects BIAS ----
gg_peak_size_bias_abs <- ggplot(df_out_inf_all_control_summ %>% 
                                  filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.000 & time <= 70 & 
                                           n_hhsize %in% c(1, 3, 5) &
                                           eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi != 1), 
                                aes(x = `Household size`, y = max_Inftot_diff_bias_abs, 
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
gg_peak_size_bias_abs

### Time of maximum infections
gg_peak_time_bias_rel <- ggplot(df_out_inf_all_control_summ %>% 
                                  filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.000 & time <= 70 & 
                                           n_hhsize %in% c(1, 3, 5) &
                                           eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi != 1) %>%
                                  # Remove Inf values
                                  mutate(max_Inftot_time_diff_bias_rel = ifelse(is.infinite(max_Inftot_time_diff_bias_rel), 
                                                                                NaN, 
                                                                                max_Inftot_time_diff_bias_rel)), 
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
# ggsave(plot = gg_peak_time_bias_rel, 
#        filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_peak_time_bias_rel.pdf", 
#        width = 12, height = 8)
# ggsave(plot = gg_peak_time_bias_rel, 
#        filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_peak_time_bias_rel.png", 
#        width = 12, height = 8)
# ggsave(plot = gg_peak_time_bias_rel, 
#        filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_peak_time_bias_rel.jpeg", 
#        width = 12, height = 8)

gg_peak_time_bias_abs <- ggplot(df_out_inf_all_control_summ %>% 
                                  filter(r_beta == 0.25 & r_tau == 0.40 & r_omega == 0.000 & time <= 70 & 
                                           n_hhsize %in% c(1, 3, 5) &
                                           eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi != 1), 
                                aes(x = `Household size`, y = max_Inftot_time_diff_bias_abs, 
                                    group = Isize, alpha = Isize, fill = Esize)) + # linetype = as.factor(eff_vax))
  geom_bar(stat = "identity", position = position_dodge(0.8)) +
  facet_grid(NPIeff ~ Esize) +
  scale_fill_discrete(l = 50) +
  scale_alpha_manual(values = alphas) +
  # xlab("Household size") +
  scale_y_continuous("Absolute bias on peak time (days)", 
                     # breaks = c(0, 250, 500, 750, 1000, 1250), 
                     labels = function(x) scales::comma(x)) +
  theme_bw(base_size = 16) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(hjust = 0, face = "bold", size = 12),
        legend.position = c(.79, .86),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank())
gg_peak_time_bias_abs

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
# ggsave(plot = gg_epidemic_size_bias_rel, 
#        filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_epidemic_size_bias_rel.pdf", 
#        width = 12, height = 8)
# ggsave(plot = gg_epidemic_size_bias_rel, 
#        filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_epidemic_size_bias_rel.png", 
#        width = 12, height = 8)
# ggsave(plot = gg_epidemic_size_bias_rel, 
#        filename = "figs/SMDM_2021/SMDM_hh_MC_SEIR_control_measures_epidemic_size_bias_rel.jpeg", 
#        width = 12, height = 8)

### Table 3: Meta regression on category-specific outcomes on Control Measures ----
df_control_measures_bias <- df_out_inf_all_control_summ %>% 
  filter(r_omega == 0.000 & 
           !((`Multicompartment structure` %in% c("E=1, I=1") & n_hhsize %in% c(1))) & 
           eff_vax %in% c(1.0) & vax_prop == 0.0, level_npi %in% c(0.4)) 

## Without interactions
fit_hh_peak_time_bias_rel <- lm(max_Inftot_time_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                  n_hhsize + 
                                  r_tau + r_beta, 
                                data = df_control_measures_bias %>% 
                                  filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel)

fit_hh_peak_size_bias_rel <- lm(max_Inftot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                                  n_hhsize + 
                                  r_tau + r_beta,
                                data = df_control_measures_bias)
summary(fit_hh_peak_size_bias_rel)

fit_hh_epidemic_size_rel <- lm(CumInfTot_diff_bias_rel ~ n_exp_states + n_inf_states + 
                             n_hhsize + r_tau + r_beta, 
                           data = df_control_measures_bias)
summary(fit_hh_epidemic_size_rel)
## With interactions
fit_hh_peak_time_bias_rel_int <- lm(max_Inftot_time_diff_bias_rel ~ n_exp_states*n_inf_states + 
                             n_exp_states*n_hhsize + 
                             n_inf_states*n_hhsize + 
                             r_tau + r_beta, 
                           data = df_control_measures_bias %>% 
                             filter(max_Inftot_time_diff_bias_rel != Inf))
summary(fit_hh_peak_time_bias_rel_int)

fit_hh_peak_size_bias_rel_int <- lm(max_Inftot_diff_bias_rel ~ n_exp_states*n_inf_states + 
                        n_exp_states*n_hhsize + 
                        n_inf_states*n_hhsize +
                        r_tau + r_beta,
                      data = df_control_measures_bias)
summary(fit_hh_peak_size_bias_rel_int)

fit_hh_epidemic_size_rel_int <- lm(CumInfTot_diff_bias_rel ~ n_exp_states*n_inf_states + 
                        n_exp_states*n_hhsize + 
                        n_inf_states*n_hhsize + 
                        r_tau + r_beta, 
                      data = df_control_measures_bias)
summary(fit_hh_epidemic_size_rel_int)

stargazer(fit_hh_peak_time_bias_rel, fit_hh_peak_time_bias_rel_int,
          fit_hh_peak_size_bias_rel, fit_hh_peak_size_bias_rel_int,
          fit_hh_epidemic_size_rel, fit_hh_epidemic_size_rel_int,
          type = "text", 
          title = "Metaregression estimates",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          covariate.labels = c("E", "I", "Household size (HH)", "$\\tau$", "$\\beta$", "E*I", "E*HH", "I*HH"), # 
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE)

stargazer(fit_hh_peak_time_bias_rel, fit_hh_peak_time_bias_rel_int,
          fit_hh_peak_size_bias_rel, fit_hh_peak_size_bias_rel_int,
          fit_hh_epidemic_size_rel, fit_hh_epidemic_size_rel_int,
          type = "latex", 
          out = "output/Table3_metaregression_bias.tex", 
          title = "Metaregression estimates",
          column.labels = c("Peak time", "Peak size", "Epidemic size"), 
          column.separate = c(2, 2, 2), 
          dep.var.caption = "", 
          dep.var.labels.include = FALSE, 
          covariate.labels = c("E", "I", "Household size (HH)", "$\\tau$", "$\\beta$", "E*I", "E*HH", "I*HH"), # 
          model.numbers = FALSE,
          omit.stat = c("all"), report = c("vc*"),
          align = TRUE)


### Appendix Figures ----
#### Appendix Figure 1 ----
max_Inftot_time_NPI20_E1_I1_hh1 <- df_out_inf_all_control_summ %>% 
  filter(r_beta == 0.25 & r_tau == 0.50 & 
           r_omega == 0.000 & time <= 70 & 
           `Multicompartment structure` %in% c("E=1, I=1") & 
           n_hhsize %in% 1 & 
           eff_vax == 1 & 
           vax_prop == 0 & 
           level_npi %in% c(0.8)) %>%
  select(max_Inftot_time, max_Inftot)
max_Inftot_time_NPI20_E1_I1_hh1$max_Inftot <- max_Inftot_time_NPI20_E1_I1_hh1$max_Inftot/10e6

max_Inftot_time_NPI20_E1_I1_hh3 <- df_out_inf_all_control_summ %>% 
  filter(r_beta == 0.25 & r_tau == 0.50 & 
           r_omega == 0.000 & time <= 70 & 
           `Multicompartment structure` %in% c("E=1, I=1") & 
           n_hhsize %in% 3 & 
           eff_vax == 1 & 
           vax_prop == 0 & 
           level_npi %in% c(0.8)) %>%
  select(max_Inftot_time, max_Inftot)
max_Inftot_time_NPI20_E1_I1_hh3$max_Inftot <- max_Inftot_time_NPI20_E1_I1_hh3$max_Inftot/10e6

max_Inftot_time_NPI20_E1_I1 <- df_out_inf_all_control_summ %>% 
  filter(r_beta == 0.25 & r_tau == 0.50 & 
           r_omega == 0.000 & time <= 70 & 
           `Multicompartment structure` %in% c("E=1, I=1", 
                                               "E=3, I=3", 
                                               "E=3, I=1", 
                                               "E=1, I=3") & 
           n_hhsize %in% c(1, 3) & 
           eff_vax == 1 & 
           vax_prop == 0 & 
           level_npi %in% c(0.8)) %>%
  mutate(x0 = c(rep(max_Inftot_time_NPI20_E1_I1_hh1$max_Inftot_time, 4),
                rep(max_Inftot_time_NPI20_E1_I1_hh3$max_Inftot_time, 4))) %>% 
  mutate(y0 = c(rep(max_Inftot_time_NPI20_E1_I1_hh1$max_Inftot, 4),
                rep(max_Inftot_time_NPI20_E1_I1_hh3$max_Inftot, 4))) %>% 
  mutate(origin = (x0 == max_Inftot_time),
         x0 = ifelse(origin, NA, x0),
         y0 = ifelse(origin, NA, y0)) %>%
  select(`Household size`,
         `Multicompartment structure`, 
         max_Inftot_time, 
         max_Inftot,
         origin,
         x0, y0)

max_Inftot_time_NPI20_E1_I1$max_Inftot <- max_Inftot_time_NPI20_E1_I1$max_Inftot/10e6

df_figA1 <- left_join(df_out_inf_all %>% 
                        filter(r_beta == 0.25 & r_tau == 0.50 & 
                                 r_omega == 0.000 & time <= 60 & 
                                 `Multicompartment structure` %in% c("E=1, I=1", 
                                                                     "E=3, I=3", 
                                                                     "E=3, I=1", 
                                                                     "E=1, I=3") & 
                                 n_hhsize %in% c(1, 3) & 
                                 eff_vax == 1 & 
                                 vax_prop == 0 & 
                                 level_npi %in% c(0.8)),
                      max_Inftot_time_NPI20_E1_I1)
gg_epidemic_curve_NPI20_E1_I1_E3_I3 <- ggplot(df_figA1, 
                                              aes(x = time, y = Inftot/10e6, color = `Household size`)) + # 
  geom_line(size = 1.1) +
  geom_segment(aes(x = max_Inftot_time_NPI20_E1_I1_hh1$max_Inftot_time, 
                   xend = max_Inftot_time_NPI20_E1_I1_hh1$max_Inftot_time, 
                   y = 0, 
                   yend = max_Inftot_time_NPI20_E1_I1_hh1$max_Inftot),
               linetype = "dashed", color = "blue") +
  geom_segment(aes(x = max_Inftot_time_NPI20_E1_I1_hh3$max_Inftot_time, 
                   xend = max_Inftot_time_NPI20_E1_I1_hh3$max_Inftot_time, 
                   y = 0, 
                   yend = max_Inftot_time_NPI20_E1_I1_hh3$max_Inftot),
               linetype = "dashed", color = "red") +
  geom_point(aes(x = max_Inftot_time_NPI20_E1_I1_hh1$max_Inftot_time, 
                 y = max_Inftot_time_NPI20_E1_I1_hh1$max_Inftot),
             color = "blue", size = 2.5) +
  geom_point(aes(x = max_Inftot_time_NPI20_E1_I1_hh3$max_Inftot_time, 
                 y = max_Inftot_time_NPI20_E1_I1_hh3$max_Inftot),
             color = "red", size = 2.5) +
  geom_segment(aes(x = x0, xend = max_Inftot_time,
                   y = y0, yend = max_Inftot,
                   color = `Household size`),
               linejoin = "mitre",
               linetype = "dashed", arrow = arrow(length = unit(0.15, "inches"),
                                                  type = "closed"), 
               show.legend = FALSE) +
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
gg_epidemic_curve_NPI20_E1_I1_E3_I3
ggsave(plot = gg_epidemic_curve_NPI20_E1_I1_E3_I3, 
       filename = "figs/Paper/FigA1_NPI20_curves.pdf", 
       width = 12, height = 8)


#### Appendix Figure 2 ----
max_Inftot_time_NPI60_E1_I1_hh1 <- df_out_inf_all_control_summ %>% 
  filter(r_beta == 0.25 & r_tau == 0.50 & 
           r_omega == 0.000 & time <= 70 & 
           `Multicompartment structure` %in% c("E=1, I=1") & 
           n_hhsize %in% 1 & 
           eff_vax == 1 & 
           vax_prop == 0 & 
           level_npi %in% c(0.4)) %>%
  select(max_Inftot_time, max_Inftot)
max_Inftot_time_NPI60_E1_I1_hh1$max_Inftot <- max_Inftot_time_NPI60_E1_I1_hh1$max_Inftot/10e6

max_Inftot_time_NPI60_E1_I1_hh3 <- df_out_inf_all_control_summ %>% 
  filter(r_beta == 0.25 & r_tau == 0.50 & 
           r_omega == 0.000 & time <= 70 & 
           `Multicompartment structure` %in% c("E=1, I=1") & 
           n_hhsize %in% 3 & 
           eff_vax == 1 & 
           vax_prop == 0 & 
           level_npi %in% c(0.4)) %>%
  select(max_Inftot_time, max_Inftot)
max_Inftot_time_NPI60_E1_I1_hh3$max_Inftot <- max_Inftot_time_NPI60_E1_I1_hh3$max_Inftot/10e6

max_Inftot_time_NPI60_E1_I1 <- df_out_inf_all_control_summ %>% 
  filter(r_beta == 0.25 & r_tau == 0.50 & 
           r_omega == 0.000 & time <= 70 & 
           `Multicompartment structure` %in% c("E=1, I=1", 
                                               "E=3, I=3", 
                                               "E=3, I=1", 
                                               "E=1, I=3") & 
           n_hhsize %in% c(1, 3) & 
           eff_vax == 1 & 
           vax_prop == 0 & 
           level_npi %in% c(0.4)) %>%
  mutate(x0 = c(rep(max_Inftot_time_NPI60_E1_I1_hh1$max_Inftot_time, 4),
                rep(max_Inftot_time_NPI60_E1_I1_hh3$max_Inftot_time, 4))) %>% 
  mutate(y0 = c(rep(max_Inftot_time_NPI60_E1_I1_hh1$max_Inftot, 4),
                rep(max_Inftot_time_NPI60_E1_I1_hh3$max_Inftot, 4))) %>% 
  mutate(origin = (x0 == max_Inftot_time),
         x0 = ifelse(origin, NA, x0),
         y0 = ifelse(origin, NA, y0)) %>%
  select(`Household size`,
         `Multicompartment structure`, 
         max_Inftot_time, 
         max_Inftot,
         origin,
         x0, y0)

max_Inftot_time_NPI60_E1_I1$max_Inftot <- max_Inftot_time_NPI60_E1_I1$max_Inftot/10e6

df_figA2 <- left_join(df_out_inf_all %>% 
                        filter(r_beta == 0.25 & r_tau == 0.50 & 
                                 r_omega == 0.000 & time <= 60 & 
                                 `Multicompartment structure` %in% c("E=1, I=1", 
                                                                     "E=3, I=3", 
                                                                     "E=3, I=1", 
                                                                     "E=1, I=3") & 
                                 n_hhsize %in% c(1, 3) & 
                                 eff_vax == 1 & 
                                 vax_prop == 0 & 
                                 level_npi %in% c(0.4)),
                      max_Inftot_time_NPI60_E1_I1)
gg_epidemic_curve_NPI60_E1_I1_E3_I3 <- ggplot(df_figA2, 
                                              aes(x = time, y = Inftot/10e6, color = `Household size`)) + # 
  geom_line(size = 1.1) +
  geom_segment(aes(x = max_Inftot_time_NPI60_E1_I1_hh1$max_Inftot_time, 
                   xend = max_Inftot_time_NPI60_E1_I1_hh1$max_Inftot_time, 
                   y = 0, 
                   yend = max_Inftot_time_NPI60_E1_I1_hh1$max_Inftot),
               linetype = "dashed", color = "blue") +
  geom_segment(aes(x = max_Inftot_time_NPI60_E1_I1_hh3$max_Inftot_time, 
                   xend = max_Inftot_time_NPI60_E1_I1_hh3$max_Inftot_time, 
                   y = 0, 
                   yend = max_Inftot_time_NPI60_E1_I1_hh3$max_Inftot),
               linetype = "dashed", color = "red") +
  geom_point(aes(x = max_Inftot_time_NPI60_E1_I1_hh1$max_Inftot_time, 
                 y = max_Inftot_time_NPI60_E1_I1_hh1$max_Inftot),
             color = "blue", size = 2.5) +
  geom_point(aes(x = max_Inftot_time_NPI60_E1_I1_hh3$max_Inftot_time, 
                 y = max_Inftot_time_NPI60_E1_I1_hh3$max_Inftot),
             color = "red", size = 2.5) +
  geom_segment(aes(x = x0, xend = max_Inftot_time,
                   y = y0, yend = max_Inftot,
                   color = `Household size`),
               linejoin = "mitre",
               linetype = "dashed", arrow = arrow(length = unit(0.15, "inches"),
                                                  type = "closed"), 
               show.legend = FALSE) +
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
gg_epidemic_curve_NPI60_E1_I1_E3_I3
ggsave(plot = gg_epidemic_curve_NPI60_E1_I1_E3_I3, 
       filename = "figs/Paper/FigA2_NPI60_curves.pdf", 
       width = 12, height = 8)