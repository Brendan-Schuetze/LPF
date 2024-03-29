##################################
#      Sim Helper Functions      #
##################################

#Only uncomment this to test using non-logistic relationship
# logistic <- function(x, a, c, d, z = NA) {
# 
#   output = pmax(0, pmin(1, if_else(x < d,
#                                    true = c,
#                                    false = (x-d) * a + c)))
# 
#   return(output)
# 
# }

# print number of parameter combinations
countParams <- function(param_list) {
  for(i in 1:length(param_list)) {
    if(i == 1) {
      poss_count <- length(unique(param_list[[i]])) * 1.0
    } else {
      poss_count <- poss_count * length(unique(param_list[[i]]))
    }
  }
  return(poss_count)
}


# Clamp between 0 and 1
clamp01 <- function(x) {
  return(pmin(0.99, pmax(0, x)))
}

# This is used to count the number of interventions applied
notEqualsZero <- function(x) {
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

logisFx <- function(x, dat_sim, summary_stat = FALSE) {
  
  logis_a <- dat_sim$Logis_a
  logis_c <- dat_sim$Logis_c
  logis_d <- dat_sim$Logis_d
  logis_z <- dat_sim$NonTreat_logis_z
  
  # If we just want to see a single curve for the means
  if(summary_stat) {
    logistic(x, 
             d = mean(logis_d), 
             z = mean(logis_z), 
             a = mean(logis_a), 
             c = mean(logis_c))
    
  } else if(length(logis_a) > 1 | 
            length(logis_d) > 1 | 
            length(logis_z) > 1 | 
            length(logis_c) > 1) {
    
    if(length(x) == sim_params$n_obs[i]) {
      xy <- list(x)
    } else {
      xy <- x
    }
    
    results_list <- lapply(xy , function(x){ 
      logistic(x,
               d = logis_d, 
               z = logis_z, 
               a = logis_a, 
               c = logis_c)})
    
    results <- as.data.frame(results_list)
    
    return(results)
    
  } else {
    logistic(x, 
             d = logis_d, 
             z = logis_z, 
             a = logis_a, 
             c = logis_c)
  }
  
}

# If performance at Time_Limit is greater than than performance criterion
# calculate time to reach performance criterion
getTimesFast_Multi <- function(dat_sim) {
  
  PC <- dat_sim$Logis_PC
  Treat <- dat_sim$Treat
  max_TL <- max(dat_sim$TL)
  
  # Generate Logis Curves
  baseline_seq <- round(seq(0, max_TL * 2, by = .10), 1)
  baseline_out <- logisFx(baseline_seq, dat_sim)
  
  times <- rep(NA, length(PC))
  
  # Go along curve and if PC is reached return time
  for(i in 1:ncol(baseline_out)) {
    times[PC <= baseline_out[, i] & is.na(times)] <- baseline_seq[i]
    
    # Return output once all times have been determined
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
