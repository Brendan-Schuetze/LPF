##############################
#       Load Libraries       #
##############################

library(MASS)
library(dplyr)
library(magrittr)
library(ggplot2)
library(psych)
library(tictoc)
library(foreach)
library(doParallel)
library(purrr)
library(ineq)
library(truncnorm)

# Import Custom Built Helper Functions
source("simHelper.R")
source("simulateIntervention.R")

# This data frame collects the summary of each simulated student population characteristics,
# and intervention's effects
sim_summary <- data.frame()

##############################
#     Simulation Settings    #
##############################

# Export data to disk for analysis?
write_data <- TRUE

# How many different combinations of parameters to test
# Simulate a little over 100,000 target because some runs will discarded

n_comb <- 10000

# How many simulated students per combination
n_obs <- 10000

# How many sims per chunk?
n_chunk <- 1000

#######################################
#   Set Values for Simulation Inputs  #
#######################################

param_list <- list(
  
  # Baseline Motivation Characteristics
  NonTreat_PC_TL_Cor = seq(0, 0.85, by = 0.01),
  NonTreat_C_Mot_Cor = c(0, 0), # Constant
  
  NonTreat_PC_Shape1 = seq(0.5, 14, by = 0.25),
  NonTreat_PC_Shape2 = seq(1, 3, by = 0.25),
  
  NonTreat_TL_Mean = seq(0.25, 8, by = 0.25),
  NonTreat_TL_sd = seq(0.05, 4, by = 0.25),
  
  NonTreat_MC_Bias = seq(-0.20, 0.20, by = 0.025),
  NonTreat_MC_sd = seq(0.02, 0.12, by = .005),
  
  NonTreat_Baseline_sd = seq(0.005, 0.10, by = 0.005),
  
  # Baseline Logistic Curve Characteristics
  NonTreat_logis_d = seq(-1, 6, by = 0.1),   # Normally distributed
  NonTreat_logis_d_sd = seq(0.01, 0.25, by = 0.05),
  
  NonTreat_logis_z = c(1, 1),  # Constant
  
  NonTreat_logis_a = seq(0.25, 1.5, by = 0.05),   # Converted to Log-Normal
  NonTreat_logis_a_sd = seq(0.05, 0.5, by = 0.05),
  
  NonTreat_logis_c_Shape1 = seq(0.25, 5, by = 0.25),   # Beta
  NonTreat_logis_c_Shape2 = seq(2, 4, by = 0.1),
  
  # Treatment effects on logistic curve
  Treat_logis_diff_a = c(rep(0, 9), seq(0, 1, by = 0.20)),
  Treat_logis_diff_c = c(rep(0, 20), seq(0, 0.5, by = 0.05)),
  Treat_logis_diff_d = c(rep(0, 80), seq(-2, 0, by = 0.05)),
  
  # Does the intervention affect motivation in terms of criterion or time?
  # We are overweighting zeroes to have more interventions where no motivation difference is applied
  Treat_Motivation_diff_PC = c(rep(0, 50), seq(0, 0.25, by = 0.01)), 
  Treat_Motivation_diff_Time = c(rep(0, 75), seq(0, 4, by = 0.1)),
  
  # Treatment metacognitive effects
  # We are overweighting zeroes to have more interventions where no metacog difference is applied
  Treat_MC_diff_Bias = c(rep(0, 15), seq(-0.10, 0.10, by = 0.025)),
  Treat_MC_diff_sd = c(rep(0, 8), seq(-0.05, 0.05, by = 0.025))
  
)

# Calculate total possible number of combinations
countParams(param_list)

# Map list to dataframe
sim_params <- map_dfr(
  param_list,
  ~ sample(.x, size = n_comb, replace = TRUE)
  )

# Add some additional metadata to the sim_params dataframe
sim_params$n_obs <- n_obs

# Count number of interventions applied
sim_params$Num_Interventions <- sim_params %>% 
  select(starts_with("Treat")) %>% 
  mutate_all(notEqualsZero) %>% 
  rowSums()

sim_params$No_Intervention_flag <- if_else(sim_params$Num_Interventions == 0,
                                           true = "No Intervention",
                                           false = "Intervention upon at least one factor")

# Filter out parameter combinations with impossible SDs
sim_params <- sim_params %>% 
  filter((NonTreat_MC_sd + Treat_MC_diff_sd) > 0) %>%
  mutate(row_n = 1:n())

##############################
#  Main Simulation Process   #
##############################

tic() # Start timer to track performance

my.cluster <- parallel::makeCluster(
  spec = detectCores() - 3, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

# Calculate sim chunks
n_sims <- nrow(sim_params)
n_loops <- ceiling(n_sims / n_chunk)

for(i in 1:n_loops) {
  
  # Determine starting and stop position
  start_n <- n_chunk * (i - 1)
  end_n <- pmin(n_sims, n_chunk * i)
  
  sim_summary_chunk <- foreach(i = start_n:end_n, .inorder = FALSE,
                         .combine = bind_rows,
                         .packages = c("dplyr", "psych", "MASS", "ggplot2", "effectsize", "data.table", "truncnorm"),
                         .errorhandling = "remove") %dopar% {
                           
                           sim_result <- simulateIntervention(sim_params, i)
                           return(sim_result)
                    
                         }
  
  sim_summary <- rbind(sim_summary, sim_summary_chunk)
  
  # Output Interim Data In Case Sim Crashes
  write.csv(x = sim_summary, file = paste("sim_summary_temp.csv", sep = ""))
  
}

toc()
stopImplicitCluster()

###################
#    Sim Export   #
###################

# Fix names
sim_summary_filtered <- sim_summary %>% tibble::remove_rownames() %>%
  mutate(Achievement_d_unstd = Treated_Perf_Mean - NonTreated_Perf_Mean) %>%
  mutate(Achievement_d_SE = (Achievement_d_CI_High - Achievement_d) / 1.96,
         Time_Spent_d_SE = (Time_Spent_d_CI_High - Time_Spent_d) / 1.96)

colnames(sim_summary_filtered) <- gsub(colnames(sim_summary_filtered), 
                              pattern = " ",
                              replacement = "_") %>%
  gsub(x = ., pattern = "^0|X0", replacement = "NonTreated") %>%
  gsub(x = ., pattern = "^1|X1", replacement = "Treated") %>%
  gsub(x = ., pattern = "\\.", replacement = "_")

# Zero out NAs because there were no applicable participants in a condition
sim_summary_filtered %<>% 
  zeroOut("NonTreated_Stopped_due_to_Anticipated_Performance") %>% 
  zeroOut("Treated_Stopped_due_to_Anticipated_Performance") %>%
  zeroOut("NonTreated_Stopped_due_to_Time") %>%
  zeroOut("Treated_Stopped_due_to_Time")

# Export
if(write_data) {
  write.csv(sim_params, "sim_params.csv")
  write.csv(sim_summary_filtered, "sim_summary.csv")
}
