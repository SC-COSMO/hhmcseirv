rm(list= ls())

GLOBAL_BATCH_FILE_PREFIX <- "./jobs_nathist"
GLOBAL_BATCH_FILE_SUFFIX <- ".txt"

GLOBAL_RSCRIPT_PATH <- "C:\\Program Files\\R\\R-3.6.3\\bin\\x64\\Rscript.exe"


GLOBAL_RUN_PATH <- "analysis/_hh_seir_dx_nathist_run_doe.R"

GLOBAL_PARAM_FILE <- "data/df_doe_mc_seirv_nathist.RDS"
df_doe_mc_seirv_runs <- readRDS(GLOBAL_PARAM_FILE)
v_pid <- unique(df_doe_mc_seirv_runs$pid)

sink(paste0(GLOBAL_BATCH_FILE_PREFIX, GLOBAL_BATCH_FILE_SUFFIX))

for(pid in v_pid) {
  cat(paste0('\"',
             GLOBAL_RSCRIPT_PATH,
             '\" \"',
             GLOBAL_RUN_PATH, 
             '\" ',
             pid))
  cat("\n")
}
sink() 
