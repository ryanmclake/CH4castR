#################################################################
# CH4cast                                                       #
# Ryan McClure                                                  #
# JAGS MODEL and FORECASTING SCRIPT                             #
#################################################################
set.seed(329)
# The JAGS MODELs
model.lm19 = ("lm.txt")
jagsscript = cat("
model {
   mu ~ dnorm(3.23,1/0.114);  # intercet
   beta ~ dnorm(0.797,1/0.000270); # cat temp paramter
   sd.pro ~ dunif(0.00001, 1000);
   tau.pro <-  pow(sd.pro, -2)

   for(i in 2:N) {
      predX[i] <- mu + C[i]*beta;
      X[i] ~ dnorm(predX[i],tau.pro); # Process variation
      Y[i] ~ dnorm(X[i], tau.obs[i]); # Observation variation
      C[i] ~ dnorm(C_mean[i],tau.pre[i]); # Covariate variation
   }
}
",
                 file = model.lm19)


model.ar19 = ("ar2_model.txt")
jagsscript = cat("
model {  
   
   #priors===================================================
   
   mu2 ~ dnorm(-0.426,1/0.0942);
   sd.pro ~ dunif(0.0001, 1000);
   tau.pro <-  pow(sd.pro, -2);
   phi ~ dnorm(0.922,1/0.00120);
   omega ~ dnorm(0.0306,1/0.000307);
   
   #Informative priors on initial conditions based on first observation
   predY[1] <- X[1];
   Y[1] ~ dnorm(X[1], tau.obs[1]);
   
   #end priors===============================================
   
   for(i in 2:N) {
      
      #process model=============================================
      
      predX[i] <- mu2 + phi*X[i-1] + omega*D[i];
      X[i] ~ dnorm(predX[i],tau.pro);
      
      #end of process model======================================
      
      #data model================================================
      
      Y[i] ~ dnorm(X[i], tau.obs[i]); # Observation variation
      
      #end of data model=========================================
   }
   IC_forecast <- X[N]}", file = model.ar19)


# # Dates to forecast in 2019 --> Based off of dates starting from when the day ebullition was measured
dates <- c(as.Date("2019-06-03"),
           as.Date("2019-06-10"),
           as.Date("2019-06-17"),
           as.Date("2019-06-24"),
           as.Date("2019-07-01"),
           as.Date("2019-07-08"),
           as.Date("2019-07-15"),
           as.Date("2019-07-22"),
           as.Date("2019-07-29"),
           as.Date("2019-08-05"),
           as.Date("2019-08-12"),
           as.Date("2019-08-19"),
           as.Date("2019-08-28"),
           as.Date("2019-09-02"),
           as.Date("2019-09-11"),
           as.Date("2019-09-20"),
           as.Date("2019-09-27"),
           as.Date("2019-10-02"),
           as.Date("2019-10-11"),
           as.Date("2019-10-16"),
           as.Date("2019-10-23"),
           as.Date("2019-10-30"),
           as.Date("2019-11-07"))

 
t1 <- proc.time()

# Sequence through the dates and the traps and execute the JAGS model

for(s in 1:length(dates)){

    # Select site and forecast date (#hashed out lines are for testing purposes)
    full_ebullition_model_alltrap_jags <- full_ebullition_model_alltrap %>% 
      filter(time <= as.Date(dates[s]))%>%
      #filter(time <= as.Date("2019-06-03"))%>%
      filter(time >= as.Date("2019-05-27"))%>%
      mutate(hobo_temp = ifelse(time<=as.Date("2019-06-17"),rnorm(length(.),mean = mean(full_ebullition_model_17[21:42,2]),sd = 1),hobo_temp),
             hobo_temp_sd = ifelse(time<=as.Date("2019-06-17"),rnorm(length(.),mean = mean(full_ebullition_model_17[21:42,3]),sd = 0.5),hobo_temp_sd))%>%
      arrange(time)
    
    
    for (i in colnames(full_ebullition_model_alltrap_jags[,c(2:10,12)])) {
      full_ebullition_model_alltrap_jags[,i] <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags[,i],option = "linear")   
    }
    
    # start the data frame used in JAGS on the first day ebullition was collected from the upstream traps
    y_nogaps <- full_ebullition_model_alltrap_jags$log_ebu_rate[!is.na(full_ebullition_model_alltrap_jags$log_ebu_rate)]
    sd_no_gaps <- full_ebullition_model_alltrap_jags$log_ebu_rate_sd[!is.na(full_ebullition_model_alltrap_jags$log_ebu_rate_sd)]
    
    # read the .nc file[h,s] from FLARE runs that correspond directly to dates ebullition from the traps were measured
    files <- list.files("./forecast_flare_driver/")                                          
    start_forecast <- max(full_ebullition_model_alltrap_jags$time)         
    start <- paste(gsub("-", "_", start_forecast))
    forecast_file <- grep(start, files, value=TRUE)
    nc <- nc_open(paste0("./forecast_flare_driver/",forecast_file))
    t <- ncvar_get(nc,'time')
    time <- as.POSIXct(t, origin = '1970-01-01 00:00.00 UTC', tz = "EST")
    time <- strftime(time, format="%Y-%m-%d")
    temp <- ncvar_get(nc,'temp')
    
    #prepare forecast temperatures from FLARE to be appended to the data frame recognized by jags
    forecast_jags <- as.data.frame(cbind(time,temp[1:18,1:210,5], temp[1:18,1:210,6], temp[1:18,1:210,7],temp[1:18,1:210,8], temp[1:18,1:210,9], temp[1:18,1:210,10]))%>%
      filter(time>=start_forecast)%>%
      filter(time<=start_forecast+16)%>%
      melt(., id = 'time')%>%
      group_by(time)%>%
      mutate(value = as.numeric(value))%>%
      summarize(water_temp_dam = mean(value, na.rm = TRUE),
                water_temp_dam_sd = sd(value, na.rm = TRUE))%>%
      arrange(time)%>%
      mutate(hobo_temp = NA)%>%
      mutate(hobo_temp_sd = NA)%>%
      mutate(log_ebu_rate = NA)%>%
      mutate(log_ebu_rate_sd = NA)%>%
      mutate(daily_temp = NA)%>%
      mutate(daily_temp_sd = NA)%>%
      mutate(daily_wind = NA)%>%
      mutate(daily_wind_sd = NA)%>%
      mutate(diff_bp = NA)%>%
      select(time, hobo_temp, hobo_temp_sd,daily_temp,daily_temp_sd,daily_wind,daily_wind_sd,diff_bp,water_temp_dam,water_temp_dam_sd,log_ebu_rate,log_ebu_rate_sd)%>%
      mutate(time = as.Date(time))
  
    
    # Prepare forecast so it can be saved as an .rds file and uploaded to github
    forecast_temp_partition <- as.data.frame(cbind(time,temp[1:18,1:210,5], temp[1:18,1:210,6], temp[1:18,1:210,7],temp[1:18,1:210,8], temp[1:18,1:210,9], temp[1:18,1:210,10]))%>%
      filter(time>=start_forecast)%>%
      filter(time<=start_forecast+16)%>%
      select(-time)%>%
      t(.)

    # Close out current NC file
    nc_close(nc)
    
    ### THIS IS THE ACTUAL JAGS MODELING PART ###
    # Append the temperature forecast from FLARE to the dataframe
    full_ebullition_model_alltrap_jags <- bind_rows(full_ebullition_model_alltrap_jags, forecast_jags[-1,])
    full_ebullition_model_alltrap_jags$hobo_temp_sd <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags$hobo_temp_sd,option = "linear")
    full_ebullition_model_alltrap_jags$log_ebu_rate_sd <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags$log_ebu_rate_sd,option = "linear")
    
    
    jags.data.lm = list(X = full_ebullition_model_alltrap_jags$hobo_temp, 
                        tau.obs = 1/(full_ebullition_model_alltrap_jags$hobo_temp_sd ^ 2),
                        N = nrow(full_ebullition_model_alltrap_jags), 
                        C_mean = full_ebullition_model_alltrap_jags$water_temp_dam,
                        tau.pre = 1/(full_ebullition_model_alltrap_jags$water_temp_dam_sd ^ 2))
    
    jags.params.lm.eval = c("sd.pro", "mu", "beta","tau.pro")
    
    
    j.lm.model   <- jags.model(file = model.lm19,
                               data = jags.data.lm,
                               n.chains = 3)
    
    eval  <- coda.samples(model = j.lm.model,
                               variable.names = jags.params.lm.eval,
                               n.iter = 10000, n.burnin = 1000)
    plot(eval)
    gelman <- gelman.diag(eval)
    gelman <- as.data.frame(bind_cols(gelman$mpsrf,as.Date(dates[s])))
    names(gelman) <- c("mpsrf","forecast_date")
    saveRDS(gelman, paste0("./forecast_output/temp_scale_gelman_diagnostics_",dates[s],".rds"))
    
    temp_out_params <- eval %>%
      spread_draws(sd.pro, mu, beta) %>%
      filter(.chain == 1) %>%
      rename(ensemble = .iteration) %>%
      ungroup()%>%
      select(sd.pro, mu, beta)
      
   
    jags.params.lm = c("sd.pro", "mu", "beta", "X", "Y", "C")
    
    jags.out   <- coda.samples(model = j.lm.model,
                               variable.names = jags.params.lm, 
                               n.iter = 10000, n.burnin = 1000)
    
    scale_out_forecast <- jags.out %>%
      spread_draws(Y[day]) %>%
      filter(.chain == 1) %>%
      rename(ensemble = .iteration) %>%
      mutate(time = full_ebullition_model_alltrap_jags$time[day]) %>%
      ungroup() %>%
      select(time, Y, ensemble)%>%
      group_by(time) %>% 
      summarise(mean = mean(Y),
                upper_90 = quantile(Y, 0.90),
                lower_90 = quantile(Y, 0.10),
                upper_80 = quantile(Y, 0.80),
                lower_80 = quantile(Y, 0.20),
                upper_70 = quantile(Y, 0.70),
                lower_70 = quantile(Y, 0.30),
                upper_60 = quantile(Y, 0.60),
                lower_60 = quantile(Y, 0.40),
                var = var(Y),
                sd = sd(Y),.groups = "drop")
    
    
    full_ebullition_model_alltrap_jags <- left_join(full_ebullition_model_alltrap_jags, scale_out_forecast[,c(1,2,12)], by = "time")%>%
        mutate(forecast_temp = ifelse(is.na(hobo_temp), mean, hobo_temp))%>%
        mutate(forecast_temp_sd = sd)
    
    
    # Generate the data frame that is recognized by the jags model
    jags.data.ar = list(X = full_ebullition_model_alltrap_jags$log_ebu_rate, 
                        tau.obs = 1/(full_ebullition_model_alltrap_jags$log_ebu_rate_sd ^ 2),
                        N = nrow(full_ebullition_model_alltrap_jags), 
                        D = full_ebullition_model_alltrap_jags$forecast_temp)
    
    # Specify the parameters of interest from the model
    jags.params.ar = c("sd.pro", "mu2", "phi", "omega", "X", "Y", "D", "IC_forecast")

    #Initialize parameters 
    nchain = 3
    chain_seeds <- c(200,800,1400)
    init <- list()
    for(i in 1:nchain){
      init[[i]] <- list(sd.pro = sd(diff(y_nogaps)),
                        .RNG.name = "base::Wichmann-Hill",
                        .RNG.seed = chain_seeds[i])
    }
    
    #Initialize JAGS model
    j.model   <- jags.model(file = model.ar19,
                             data = jags.data.ar,
                             inits = init,
                             n.chains = 3)
    
    #Run JAGS model and sample from the posteriors
    eval_ebu  <- coda.samples(model = j.model,
                               variable.names = c("sd.pro", "mu2", "phi", "omega"),
                               n.iter = 10000, n.burnin = 3000)

    plot(eval_ebu)
    gelman_ebu <- gelman.diag(eval_ebu)
    gelman_ebu <- as.data.frame(bind_cols(gelman_ebu$mpsrf,as.Date(dates[s])))
    names(gelman_ebu) <- c("mpsrf","forecast_date")
    
    saveRDS(gelman_ebu, paste0("./forecast_output/ebu_model_gelman_diagnostics_",dates[s],".rds"))
    
    #Run JAGS model and sample from the posteriors
    jags.out   <- coda.samples(model = j.model,
                               variable.names = jags.params.ar, 
                               n.iter = 10000, n.burnin = 1000)
    
    ##GENEARTE THE CREDIBLE INTERVAL
    ebu_out_forecast <- jags.out %>%
      spread_draws(Y[day]) %>%
      filter(.chain == 1) %>%
      rename(ensemble = .iteration) %>%
      mutate(time = full_ebullition_model_alltrap_jags$time[day]) %>%
      ungroup() %>%
      select(time, Y, ensemble)%>%
      group_by(time) %>% 
      summarise(mean = mean(Y),
                upper_90 = quantile(Y, 0.90),
                lower_90 = quantile(Y, 0.10),
                upper_80 = quantile(Y, 0.80),
                lower_80 = quantile(Y, 0.20),
                upper_70 = quantile(Y, 0.70),
                lower_70 = quantile(Y, 0.30),
                upper_60 = quantile(Y, 0.60),
                lower_60 = quantile(Y, 0.40),
                var = var(Y),
                sd = sd(Y),.groups = "drop")
    
    # # Save just forecast from the date of the observation + 10 days into the future
    forecast_saved_ebu <- ebu_out_forecast %>%
      filter(time >= start_forecast)%>%
      mutate(forecast_date = start_forecast)
    saveRDS(forecast_saved_ebu, paste0("./forecast_output/ebullition_predictive_interval_forecast_alltrap_",dates[s],".rds"))

    # Extract the parameter esimates from the ebullition model run for dates[s] and traps[h]
    ########################################################################################
    ebu_out_parms <- jags.out %>%
      spread_draws(sd.pro, mu2, phi, omega, IC_forecast) %>%
      filter(.chain == 1) %>%
      rename(ensemble = .iteration) %>%
      mutate(forecast_date = start_forecast)%>%
      ungroup()%>%
      select(forecast_date, sd.pro, mu2, phi, omega, IC_forecast)
    saveRDS(ebu_out_parms, paste0("./forecast_output/ebullition_parameters_",dates[s],".rds"))
    #######################################################################################

    
    #UNCERTAINTY PARTITIONING OF DATA
    ########################################################################################

    Nmc <- 1000         ## set number of Monte Carlo draws
    future <- seq(from = start_forecast, to = start_forecast+lubridate::days(16), by = "days")

    
    forecast_function <- function(IC, mu2, phi, omega, FLARE, Q = 0, n = Nmc){
      N <- matrix(NA, n, 17)  ## storage
      Npre <- IC
      for(t in 1:ncol(N)){
        est <- mu2 + phi * Npre + omega * FLARE + Q
        N[,t] <- rnorm(n, est, Q)                                    ## predict next step
        Npre <- N[, t]                                               ## update IC
      }
      return(N)
    }
    
    parms <- sample_n(ebu_out_parms, Nmc, replace=TRUE)
    FLARE = scale_out_forecast[,c(1,2,12)]%>%filter(time>=start_forecast)
    drow <- sample.int(nrow(forecast_temp_partition), Nmc, replace = TRUE)
    
    # initial condition uncertainty
    IC_uncertainty <- forecast_function(IC = parms$IC_forecast,  ## sample IC
                             FLARE = mean(FLARE$mean),
                             mu2 = mean(parms$mu2),
                             phi = mean(parms$phi),
                             omega = mean(parms$omega),
                             Q = 0,
                             n = Nmc)
    
    IC_uncertainty <- t(IC_uncertainty)
    IC_uncertainty <- cbind.data.frame(future, IC_uncertainty)%>%
      rename(time = future)
    
    IC_uncertainty <- melt(IC_uncertainty, id = "time")
    
    IC_partition <- IC_uncertainty %>%
      group_by(time) %>%
      summarise(mean = mean(value),
                max = max(value),
                min = min(value),
                upper = quantile(value, 0.95),
                lower = quantile(value, 0.05),
                var = var(value),
                sd = sd(value))%>%
      mutate(forecast_date = start_forecast)
    saveRDS(IC_partition, paste0("./forecast_output/initial_condition_",dates[s],".rds"))
    
    # parameter uncertainty
    param_uncertainty <- forecast_function(IC = mean(parms$IC_forecast),  ## sample IC
                                        FLARE = mean(FLARE$mean),
                                        mu2 = parms$mu2,
                                        phi = parms$phi,
                                        omega = parms$omega,
                                        Q = 0,
                                        n = Nmc)
    
    param_uncertainty <- t(param_uncertainty)
    param_uncertainty <- cbind.data.frame(future, param_uncertainty)%>%
      rename(time = future)
    
    param_uncertainty <- melt(param_uncertainty, id = "time")
    
    param_partition <- param_uncertainty %>%
      group_by(time) %>%
      summarise(mean = mean(value),
                max = max(value),
                min = min(value),
                upper = quantile(value, 0.95),
                lower = quantile(value, 0.05),
                var = var(value),
                sd = sd(value))%>%
      mutate(forecast_date = start_forecast)
    saveRDS(param_partition, paste0("./forecast_output/model_parameter_",dates[s],".rds"))
    
    
    # driver data uncertainty
    driver_uncertainty <- forecast_function(IC = mean(parms$IC_forecast),  ## sample IC
                                           FLARE = as.numeric(forecast_temp_partition[drow,]),
                                           mu2 = mean(parms$mu2),
                                           phi = mean(parms$phi),
                                           omega = mean(parms$omega),
                                           Q = 0,
                                           n = Nmc)
    
    driver_uncertainty <- t(driver_uncertainty)
    driver_uncertainty <- cbind.data.frame(future, driver_uncertainty)%>%
      rename(time = future)
    
    driver_uncertainty <- melt(driver_uncertainty, id = "time")
    
    drive_partition <- driver_uncertainty %>%
      group_by(time) %>%
      summarise(mean = mean(value),
                max = max(value),
                min = min(value),
                upper = quantile(value, 0.95),
                lower = quantile(value, 0.05),
                var = var(value),
                sd = sd(value))%>%
      mutate(forecast_date = start_forecast)
    saveRDS(drive_partition, paste0("./forecast_output/model_driver_",dates[s],".rds"))
    
    # process uncertainty
    process_uncertainty <- forecast_function(IC = mean(parms$IC_forecast),  ## sample IC
                                            FLARE = mean(FLARE$mean),
                                            mu2 = mean(parms$mu2),
                                            phi = mean(parms$phi),
                                            omega = mean(parms$omega),
                                            Q = sqrt(parms$sd.pro),
                                            n = Nmc)
    
    process_uncertainty <- t(process_uncertainty)
    process_uncertainty <- cbind.data.frame(future, process_uncertainty)%>%
      rename(time = future)
    
    process_uncertainty <- melt(process_uncertainty, id = "time")
    
    process_partition <- process_uncertainty %>%
      group_by(time) %>%
      summarise(mean = mean(value),
                max = max(value),
                min = min(value),
                upper = quantile(value, 0.95),
                lower = quantile(value, 0.05),
                var = var(value),
                sd = sd(value))%>%
      mutate(forecast_date = start_forecast)
    saveRDS(process_partition, paste0("./forecast_output/model_process_",dates[s],".rds"))
}

dopar_loop <- proc.time()-t1
print(dopar_loop[3]/60)
