simulateIntervention <- function(sim_params, i, return_dat_sim = FALSE) {
  
  ###############################
  # Set parameters to variables #
  ###############################
  
  NonTreat_PC_TL_Cor <- sim_params$NonTreat_PC_TL_Cor[i]
  NonTreat_C_Mot_Cor <- sim_params$NonTreat_C_Mot_Cor[i]
  
  # Generate Time Limit and Performance Cutoff
  # These are simulated at varying levels of correlation
  mu <- c(0, 0, 0)
  Sigma <- matrix(c(1, NonTreat_PC_TL_Cor, NonTreat_C_Mot_Cor,
                    NonTreat_PC_TL_Cor, 1, NonTreat_C_Mot_Cor, 
                    NonTreat_C_Mot_Cor, NonTreat_C_Mot_Cor, 1), 
                  nrow = 3, ncol = 3)  
  
  dat_sim <- data.frame(mvrnorm(n = sim_params$n_obs[i], mu = mu, Sigma = Sigma))
  colnames(dat_sim) <- c("Performance_Criterion_Raw", "Time_Limit_Raw", "Logis_C_Raw")
  
  NonTreat_PC_Shape1 <- sim_params$NonTreat_PC_Shape1[i]
  NonTreat_PC_Shape2 <- sim_params$NonTreat_PC_Shape2[i]
  NonTreat_TL_Mean <- sim_params$NonTreat_TL_Mean[i]
  NonTreat_TL_sd <- sim_params$NonTreat_TL_sd[i]
  
  NonTreat_MC_Bias <- sim_params$NonTreat_MC_Bias[i]
  dat_sim$NonTreat_MC_sd <- sim_params$NonTreat_MC_sd[i]
  NonTreat_Baseline_sd <- sim_params$NonTreat_Baseline_sd[i]
  
  # Motivation Intervention Parameters
  Treat_Motivation_diff_Time <- sim_params$Treat_Motivation_diff_Time[i]
  Treat_Motivation_diff_PC <- sim_params$Treat_Motivation_diff_PC[i]
  
  # Metacognition Intervention Parameters
  Treat_MC_diff_Bias <- sim_params$Treat_MC_diff_Bias[i]
  Treat_MC_diff_sd <- sim_params$Treat_MC_diff_sd[i]
  
  # Learning Curve Parameters
  dat_sim$Treat_logis_diff_a <- sim_params$Treat_logis_diff_a[i]
  dat_sim$Treat_logis_diff_c <- sim_params$Treat_logis_diff_c[i]
  dat_sim$Treat_logis_diff_d <- sim_params$Treat_logis_diff_d[i]
  
  # Generate Logis d from normal distribution
  dat_sim$NonTreat_logis_d <- rnorm(n = sim_params$n_obs[i],
                                    sim_params$NonTreat_logis_d[i],
                                    sim_params$NonTreat_logis_d_sd[i])
  
  # Logis z is a constant
  dat_sim$NonTreat_logis_z <- sim_params$NonTreat_logis_z[i]
  
  #NonTreat_logis_a <- sim_params$NonTreat_logis_a[i]
  dat_sim$NonTreat_logis_a <- rlnorm_norm(n = sim_params$n_obs[i], 
                                          m = sim_params$NonTreat_logis_a[i], 
                                          s = sim_params$NonTreat_logis_a_sd[i])
  
  dat_sim$Treat_logis_a <- dat_sim$NonTreat_logis_a + dat_sim$Treat_logis_diff_a
  dat_sim$Treat_logis_d <- dat_sim$NonTreat_logis_d + dat_sim$Treat_logis_diff_d
 
  # Assign Treated (Learning Curve)
  dat_sim$Treat <- sample(x = c(0, 1), size = sim_params$n_obs[i], replace = TRUE)
  
  # Assign Treated (Metacog)
  dat_sim$MC_Treat_Bias <- dat_sim$Treat * Treat_MC_diff_Bias
  dat_sim$MC_Treat_SD <- dat_sim$Treat * Treat_MC_diff_sd
  
  # Assign Treated (Motivation Intervention)
  dat_sim$Mot_Treat_Time <- dat_sim$Treat * Treat_Motivation_diff_Time
  dat_sim$Mot_Treat_PC <- dat_sim$Treat * Treat_Motivation_diff_PC
  
  ##############################
  # Inverse Transform Sampling #
  ##############################
  
  # Performance Criterion is a beta distribution
  dat_sim$PC <- qbeta(pnorm(dat_sim$Performance_Criterion_Raw, mean = 0, sd = 1), 
                      shape1 = NonTreat_PC_Shape1, 
                      shape2 = NonTreat_PC_Shape2)
  
  dat_sim$PC <- pmax(pmin(1, dat_sim$PC + dat_sim$Mot_Treat_PC), 0)
  
  # Time Limit is a Chi-Squared distribution
  # Here we are converting correlated normal dist to chi-sq
  dat_sim$TL <- qtruncnorm(pnorm(dat_sim$Time_Limit_Raw, mean = 0, sd = 1), 
                       m = NonTreat_TL_Mean + dat_sim$Mot_Treat_Time, 
                       sd = NonTreat_TL_sd, a = 0)
  
  dat_sim$NonTreat_logis_c <- qbeta(pnorm(dat_sim$Logis_C_Raw, mean = 0, sd = 1),
                                    sim_params$NonTreat_logis_c_Shape1[i], 
                                    sim_params$NonTreat_logis_c_Shape2[i])
  
  dat_sim$Treat_logis_c <- pmax(0, pmin(0.999, dat_sim$NonTreat_logis_c + dat_sim$Treat_logis_diff_c))

  #########################
  #   Adjust for Metacog  #
  #########################
  
  # Adjust Performance Criterion for Metacognitive Limits
  dat_sim$Metacog_Bias <- rnorm(n = nrow(dat_sim), 
                                mean = NonTreat_MC_Bias + dat_sim$MC_Treat_Bias, 
                                sd = dat_sim$NonTreat_MC_sd + dat_sim$MC_Treat_SD)
  
  # When metacog bias is positive (over confidence), studying stops earlier
  # when negative (underconfidence) studying stops later
  dat_sim$PC_Meta <- pmax(0, pmin(dat_sim$PC - dat_sim$Metacog_Bias, 0.999)) # MAY WANT TO SWITCH THIS ERROR TERM TO THE LOGIS TL #
  
  #############################
  # Calculate prior knowledge #
  #############################
  
  dat_sim$Prior_Knowledge <- if_else(dat_sim$Treat == 1,
                                     true = treatFx(0, dat_sim)[[1]],
                                     false = nontreatFx(0, dat_sim)[[1]])
  
  # Add the same noise to prior knowledge
  #### Could we just capture the prior knowledge test at t/2 or something?
  dat_sim$Prior_Knowledge_Meta <- pmin(dat_sim$PC_Meta, dat_sim$Prior_Knowledge)
  dat_sim$Prior_Knowledge_Bounded <- pmin(1, 
                                          pmax(0, 
                                               dat_sim$Prior_Knowledge_Meta + 
                                                    rnorm(n = sim_params$n_obs[i], 
                                                          mean = 0, sd = NonTreat_Baseline_sd)))
  
  #########################
  # Calculate Performance #
  #########################
  
  # Calculate performance if they simply studied for entire TL
  dat_sim$Logis_TL <- if_else(dat_sim$Treat == 1,
                              true = treatFx(dat_sim$TL, dat_sim)[[1]],
                              false = nontreatFx(dat_sim$TL, dat_sim)[[1]])
  
  # Would they have stopped studying earlier?
  dat_sim$Logis_PC <- if_else(dat_sim$Logis_TL < dat_sim$PC_Meta,
                              true = dat_sim$Logis_TL,
                              false = dat_sim$PC_Meta)
  
  # If Performance at Time_Limit is less than Performance criterion
  # then return performance at time limit, else calculate time to
  # reach performance criterion
  dat_sim$Time_Spent <- if_else(dat_sim$Logis_TL < dat_sim$PC_Meta,
                                true = dat_sim$TL,
                                false = getTimesFast_Multi(dat_sim))
  
  # Determine if each student stopped due to time or performance
  dat_sim$PC_TL_Stop <- if_else(dat_sim$Logis_TL < dat_sim$PC_Meta,
                                true = "Stopped due to Time",
                                false = "Stopped due to Anticipated Performance")
  
  # Make sure adjusted performance outcome is bounded between zero and one
  # Also add in baseline noise as not to have any crazy high Cohen's d
  dat_sim$Logis_Bounded <- pmin(1, pmax(0, 
                                        dat_sim$Logis_PC + rnorm(n = sim_params$n_obs[i], 
                                                                 mean = 0, sd = NonTreat_Baseline_sd)))
  
  ############################
  # Calculate Additional DVs #
  ############################
  
  # Calculate Efficiency
  # dat_sim$Efficiency <- (dat_sim$Logis_Bounded / dat_sim$Time_Spent)
  
  # Calculate Absolute and Netgain
  dat_sim$AG <- dat_sim$Logis_Bounded - dat_sim$Prior_Knowledge_Bounded
  dat_sim$NG <- (dat_sim$AG) / (1.01 - dat_sim$Prior_Knowledge_Bounded)
  
  ##########################
  # Sim Summary Statistics #
  ##########################
  
  # Generate Diagnostic Plots
  x <- seq(0, 12, by = 0.25)
  T_out <- treatFx(x, dat_sim, summary_stat = TRUE)
  NT_out <- nontreatFx(x, dat_sim, summary_stat = TRUE)
  
  # Flag Certain logistic combinations for potential pruning
  T_NT_Startpoint_Switched <- as.numeric(T_out[1] < NT_out[1])
  T_NT_Switch <- as.numeric(any(T_out < NT_out))
  
  # Export Diagnostic Graphs (treatment in red, control in black)
  if(i %% 500 == 0) {
    PC_Meta_Treat <- mean(dat_sim$PC_Meta[dat_sim$Treat == 1])
    PC_Meta_NonTreat <- mean(dat_sim$PC_Meta[dat_sim$Treat == 1])
    
    png(filename = paste0("Graphs/sim_",i,".png"))
    
    # Plot Learning Curves
    plot(type = "l", x = x, y = NT_out, xlim = c(0, 15), ylim = c(0,1),
         xlab = "Time Studying", ylab = "Performance")
    lines(x = x, y = T_out, col = "red")
    abline(v = NonTreat_TL_Mean, lty = "dotted")
    abline(v = NonTreat_TL_Mean + Treat_Motivation_diff_Time , lty = "dotted", col = "red")
    abline(h = PC_Meta_NonTreat, lty = "dotted", col = "black")
    abline(h = PC_Meta_Treat, lty = "dotted", col = "red")
    dev.off()
    
    # Export Current Simulation Status
    file_conn <- file("sim_status.txt")
    writeLines(paste("Current Row", i), file_conn)
    close(file_conn)
    
  }
  
  dat_sim$Treat <- relevel(as.factor(dat_sim$Treat), ref = "1")
  
  # Calculate Effect Sizes
  Logis_TL_d <- effectsize::cohens_d(data = dat_sim, Logis_TL ~ Treat)
  Logis_Bounded_d <- effectsize::cohens_d(data = dat_sim, Logis_Bounded ~ Treat)
  Time_Spent_d <- effectsize::cohens_d(data = dat_sim, Time_Spent ~ Treat)
  
  # | abs(Time_Spent_d$cohen.d[2] > 2) | abs(Logis_Bounded_d$cohen.d[2]) > 2
  if(i %% 500 == 0) { 
    dens_plot <- ggplot(data = dat_sim) + 
      geom_density(aes(x = Logis_Bounded, fill = as.factor(Treat)), alpha = 0.45) +
      theme_classic()
    ggsave(dens_plot, filename = paste0("Graphs/density_",i,".png"))
  }
  
  # Treated_Gini <- try(Gini(dat_sim$Logis_Bounded[dat_sim$Treat == 1]), silent = TRUE)
  # NonTreated_Gini <- try(Gini(dat_sim$Logis_Bounded[dat_sim$Treat == 0]), silent = TRUE)
  # 
  # if(class(Treated_Gini) == "try-error") {
  #   Treated_Gini <- NA
  # }
  # 
  # if(class(NonTreated_Gini) == "try-error") {
  #   NonTreated_Gini <- NA
  # }
  
  sim_results <- data.frame(
    
    # Simulation Information
    row_n = i,
    n_interventions = sim_params$Num_Interventions[i],
    no_intervention_flag = sim_params$No_Intervention_flag[i],
    n_obs = sim_params$n_obs[i],
    T_NT_Startpoint_Switched,
    T_NT_Switch,
    
    # Non-Treatment Population Parameters
    NonTreat_Baseline_sd,
    NonTreat_C_Mot_Cor,
    NonTreat_PC_TL_Cor,
    NonTreat_PC_Shape1,
    NonTreat_PC_Shape2,
    NonTreated_PC_Mean = NonTreat_PC_Shape1 / (NonTreat_PC_Shape1 + NonTreat_PC_Shape2),
    NonTreat_TL_Mean = NonTreat_TL_Mean,
    NonTreat_TL_sd = NonTreat_TL_sd,
    NonTreat_MC_Bias = sim_params$NonTreat_MC_Bias[i],
    NonTreat_MC_sd = sim_params$NonTreat_MC_sd[i],
    
    # Untreated Learning Curve Parameters (is there a better way to std.?)
    NonTreat_logis_z = sim_params$NonTreat_logis_z[i],
    
    NonTreat_logis_a = sim_params$NonTreat_logis_a[i],
    NonTreat_logis_a_sd = sim_params$NonTreat_logis_a_sd[i],
    
    NonTreat_logis_d = sim_params$NonTreat_logis_d[i],
    NonTreat_logis_d_sd = sim_params$NonTreat_logis_d_sd[i],
    
    NonTreat_logis_c_Shape1 = sim_params$NonTreat_logis_c_Shape1[i],
    NonTreat_logis_c_Shape2 = sim_params$NonTreat_logis_c_Shape2[i],
    NonTreated_logis_c_Mean = mean(dat_sim$NonTreat_logis_c),
    NonTreated_logis_c_SD = sd(dat_sim$NonTreat_logis_c),
    
    # Fundamental Treatment Effects
    Treat_logis_diff_a = sim_params$Treat_logis_diff_a[i],
    Treat_logis_diff_c = sim_params$Treat_logis_diff_c[i],
    Treat_logis_diff_d = sim_params$Treat_logis_diff_d[i],
    
    Treat_MC_diff_Bias,
    Treat_MC_diff_sd,
    Treat_Motivation_diff_Time,
    Treat_Motivation_diff_PC,
    
    # Treated Learning Curve Parameters
    Treated_logis_d = mean(dat_sim$Treat_logis_d),
    Treated_logis_a = mean(dat_sim$Treat_logis_a),
    Treated_logis_c = mean(dat_sim$Treat_logis_c),
    
    # Intervention Descriptive Statistics
    # Note that we ad "ed" to these column names to indicate that these 
    # are computed after simulation run
    Treated_Perf_Mean = mean(dat_sim$Logis_Bounded[dat_sim$Treat == 1]),
    Treated_Perf_SD = sd(dat_sim$Logis_Bounded[dat_sim$Treat == 1]),
    Treated_Prior_Knowledge_Mean = mean(dat_sim$Prior_Knowledge_Bounded[dat_sim$Treat == 1]),
    #Treated_Gini = Treated_Gini,
    
    NonTreated_Perf_Mean = mean(dat_sim$Logis_Bounded[dat_sim$Treat == 0]),
    NonTreated_Perf_SD = sd(dat_sim$Logis_Bounded[dat_sim$Treat == 0]),
    NonTreated_Prior_Knowledge_Mean = mean(dat_sim$Prior_Knowledge_Bounded[dat_sim$Treat == 0]),
    #NonTreated_Gini = NonTreated_Gini,
    
    Treated_Time_Mean = mean(dat_sim$Time_Spent[dat_sim$Treat == 1]),
    Treated_Time_SD = sd(dat_sim$Time_Spent[dat_sim$Treat == 1]),
    
    NonTreated_Time_Mean = mean(dat_sim$Time_Spent[dat_sim$Treat == 0]),
    NonTreated_Time_SD = sd(dat_sim$Time_Spent[dat_sim$Treat == 0]),
  
    # Cohen's d Outcome Measures
    Achievement_Unbounded_d_CI_Low = Logis_TL_d$CI_low,
    Achievement_Unbounded_d = Logis_TL_d$Cohens_d,
    Achievement_Unbounded_d_CI_High = Logis_TL_d$CI_high,
    
    Achievement_d_CI_Low = Logis_Bounded_d$CI_low,
    Achievement_d = Logis_Bounded_d$Cohens_d,
    Achievement_d_CI_High = Logis_Bounded_d$CI_high,
    
    Time_Spent_d_CI_Low = Time_Spent_d$CI_low,
    Time_Spent_d = Time_Spent_d$Cohens_d,
    Time_Spent_d_CI_High = Time_Spent_d$CI_high,
 
    # Simonsmeier Analyses (only using control condition data)
    Cor_PK = cor(dat_sim[dat_sim$Treat == 0,]$Prior_Knowledge_Bounded,
                 dat_sim[dat_sim$Treat == 0,]$Logis_Bounded, use = "complete.obs"),
    
    Cor_AG = cor(dat_sim[dat_sim$Treat == 0,]$Prior_Knowledge_Bounded,
                 dat_sim[dat_sim$Treat == 0,]$AG, use = "complete.obs"),
    
    Cor_NG = cor(dat_sim[dat_sim$Treat == 0,]$Prior_Knowledge_Bounded,
                 dat_sim[dat_sim$Treat == 0,]$NG,
                 use = "complete.obs")
  )
  
  ctab <- t(data.frame(table(paste(dat_sim$Treat, dat_sim$PC_TL_Stop))))
  colnames(ctab) <- ctab[1,]
  ctab <- t(data.frame(ctab[2,]))
  sim_results <- cbind(sim_results, ctab)
  
  if(return_dat_sim) {
    return(dat_sim)
  } else {
    return(sim_results)
  }
}
