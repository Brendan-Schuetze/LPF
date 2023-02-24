# Sim Helper Functions

# !! This is not actually a logistic function !!
# Only uncomment this to test using non-logistic relationship
# logistic <- function(x, a, c, d, z = NA) {
#   
#   output = pmax(0, pmin(1, if_else(x < d, 
#                                    true = c, 
#                                    false = (x-d) * a + c)))
#   
#   return(output)
#   
# }

# This is used to count the number of interventions applied
not_equals_zero <- function(x) {
  return(as.numeric(x != 0))
}

zeroOut <- function(df, x) {
  df[[x]] <- if_else(is.na(df[[x]]), true = 0, false = as.numeric(df[[x]]))
  
  return(df)
}

# Log Normal Functions that deal with non-logged means and SDs
rlnorm_norm <- function(n, m, s) {
  location <- log(m^2 / sqrt(s^2 + m^2))
  shape <- sqrt(log(1 + (s^2 / m^2)))
  return(rlnorm(n = n, meanlog = location, sdlog = shape))
}

qlnorm_norm <- function(x, m, s) {
  location <- log(m^2 / sqrt(s^2 + m^2))
  shape <- sqrt(log(1 + (s^2 / m^2)))
  return(qlnorm(p = x, meanlog = location, sdlog = shape))
}

treatFx <- function(x, dat_sim, summary_stat = FALSE) {
  
  Treat_logis_a <- dat_sim$Treat_logis_a
  Treat_logis_c <- dat_sim$Treat_logis_c
  Treat_logis_d <- dat_sim$Treat_logis_d
  Treat_logis_z <- dat_sim$NonTreat_logis_z
  
  # If we just want to see a single curve for the means
  if(summary_stat) {
    logistic(x, 
             d = mean(Treat_logis_d), 
             z = mean(Treat_logis_z), 
             a = mean(Treat_logis_a), 
             c = mean(Treat_logis_c))
  
    } else if(length(Treat_logis_a) > 1 | 
              length(Treat_logis_d) > 1 | 
              length(Treat_logis_z) > 1 | 
              length(Treat_logis_c) > 1) {
    
      if(length(x) == sim_params$n_obs[i]) {
        xy <- list(x)
      } else {
        xy <- x
      }
        
      results_list <- lapply(xy , function(x){ 
        logistic(x,
                 d = Treat_logis_d, 
                 z = Treat_logis_z, 
                 a = Treat_logis_a, 
                 c = Treat_logis_c)})
      
      results <- as.data.frame(results_list)
      
      return(results)
    
  } else {
    logistic(x, 
             d = Treat_logis_d, 
             z = Treat_logis_z, 
             a = Treat_logis_a, 
             c = Treat_logis_c)
  }
  
}

nontreatFx <- function(x, dat_sim, summary_stat = FALSE) {
  
  NonTreat_logis_a <- dat_sim$NonTreat_logis_a
  NonTreat_logis_c <- dat_sim$NonTreat_logis_c
  NonTreat_logis_d <- dat_sim$NonTreat_logis_d
  NonTreat_logis_z <- dat_sim$NonTreat_logis_z
  
  # If we just want to see a single curve for the means
  if(summary_stat) {
    logistic(x, 
             d = mean(NonTreat_logis_d), 
             z = mean(NonTreat_logis_z), 
             a = mean(NonTreat_logis_a), 
             c = mean(NonTreat_logis_c))
  
    } else if(length(NonTreat_logis_d) >1  | 
              length(NonTreat_logis_a) > 1 | 
              length(NonTreat_logis_z) > 1 | 
              length(NonTreat_logis_c) > 1) {

    if(length(x) == sim_params$n_obs[i]) {
      xy <- list(x)
    } else {
      xy <- x
    }
    
    results_list <- lapply(xy , function(x){ 
                           logistic(x,
                           d = NonTreat_logis_d, 
                           z = NonTreat_logis_z, 
                           a = NonTreat_logis_a, 
                           c = NonTreat_logis_c)})
    
    results <- as.data.frame(results_list)
    
    return(results)
    
  } else {
    logistic(x, 
           d = NonTreat_logis_d, 
           z = NonTreat_logis_z, 
           a = NonTreat_logis_a, 
           c = NonTreat_logis_c)
  }
}


# If performance at Time_Limit is greater than than performance criterion
# calculate time to reach performance criterion
getTimesFast_Multi <- function(dat_sim) {
  
  PC <- dat_sim$PC_Meta
  Treat <- dat_sim$Treat
  max_TL <- max(dat_sim$TL)
  
  # Generate Treatment and Non-Treatment Curves
  baseline_seq <- round(seq(0, max_TL * 2, by = .10), 1)
  baseline_treat <- treatFx(baseline_seq, dat_sim)
  baseline_nontreat <- nontreatFx(baseline_seq, dat_sim)
  
  times <- rep(NA, length(PC))
  
  # Go along curve and if PC is reached return time
  for(i in 1:ncol(baseline_treat)) {
    times[PC <= baseline_treat[, i] & is.na(times) & Treat == 1] <- baseline_seq[i]
    times[PC <= baseline_nontreat[, i] & is.na(times) & Treat == 0] <- baseline_seq[i]
    
    if(!any(is.na(times))) {
      return(times)
    }
  }
  
  return(times)
}

simTester <- function(sim_params, n_iter) {
  mean_ds <- c()  
  sd_ds <- c()

  for(i in 1:nrow(sim_params)) {
    ds <- c()
    for(j in 1:n_iter) {
      x2 <- simulateIntervention(sim_params, i)
      ds <- c(ds, x2$Achievement_d)
    }
    mean_ds <- c(mean_ds, mean(ds))
    sd_ds <- c(sd_ds, sd(ds))
  }
  return(data.frame(row_n = sim_params$row_n, mean_ds = mean_ds, sd_ds = sd_ds))
}
