library(dplyr)
library(ggplot2)
source("R/02_decision_model_functions.R")


GLOBAL_PARAM_FILE <- "data/df_doe_mc_seirv_nathist.RDS"
df_doe_mc_seirv_runs <- readRDS(GLOBAL_PARAM_FILE)
v_pid <- unique(df_doe_mc_seirv_runs$pid)

list_l_out <- list()
v_failed_pid <- c()
for(i in 1:length(v_pid)){ # n_pid <- v_pid[1]
  tryCatch({
    load(file = paste0("output/output_doe_mc_seirv_", v_pid[i],".RData")) 
    list_l_out[[i]] <- l_out
    names(list_l_out)[i] <- v_pid[i]
  }, error = function(cond) {
    print(paste(Sys.time(),": +++ ERROR in pid", v_pid[i], cond, "\n"))
    v_failed_pid <<- c(v_failed_pid, v_pid[i])
  })
}

v_good_pid <- as.numeric(names(list_l_out))[!is.na(as.numeric(names(list_l_out)))]
# v_failed_pid <- v_pid[is.na(as.numeric(names(list_l_out)))]

save(list_l_out, v_good_pid,
     file = "output/output_doe_mc_seirv_all_nathist.RData")

#### Analyze output from DOE Natural History ####
load(file = "output/output_doe_mc_seirv_all_nathist.RData")

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
df_out_inf_all <- df_out_inf_all %>%
  group_by(pid) %>%
  mutate(max_Inftot = max(Inftot),
         max_Inftot_time = time[which.max(Inftot)])

## Analyze epidemic outputs
# Find the t at which I(t) is at its max
# Find times at which I(t) is at a percentage of its max
# Ranges of times at which I starts to rise: compute doubling times
