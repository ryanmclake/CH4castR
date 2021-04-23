
model.ar2 = ("ar2_model.txt")
jagsscript = cat("
model {  
   
   #priors===================================================
   
   mu2 ~ dnorm(0,5);
   sd.pro ~ dunif(0.0001, 1000);
   tau.pro <-  pow(sd.pro, -2);
   phi ~ dnorm(0,5);
   omega ~ dnorm(0,5);
   delta ~ dnorm(0,5)
   
   #Informative priors on initial conditions based on first observation
   predY[1] <- X[1];
   Y[1] ~ dnorm(X[1], tau.obs[1]);
   
   #end priors===============================================
   
   for(i in 2:N) {
      
      #process model=============================================
      
      predX[i] <- mu2 + phi*X[i-1]+omega*D[i]+delta*B[i];
      X[i] ~ dnorm(predX[i],tau.pro);
      
      #end of process model======================================
      
      #data model================================================
      
      Y[i] ~ dnorm(X[i], tau.obs[i]); # Observation variation
      D[i] ~ dnorm(D_mean[i],tau.pre[i]); # temp Covariate variation
      
      #end of data model=========================================
   }
   IC_forecast <- X[N]
}
", file = model.ar2)

# # Dates to forecast in 2019 --> Based off of dates starting from when the day ebullition was measured
#dates <- seq(from = as.Date("2019-05-27"), to = as.Date("2019-11-07"), by = "days")

# # Dates to forecast in 2019 --> Based off of dates starting from when the day ebullition was measured
dates <- c(as.Date("2019-05-27"),as.Date("2019-06-03"),as.Date("2019-06-10"),
           as.Date("2019-06-17"),as.Date("2019-06-24"),as.Date("2019-07-01"),
           as.Date("2019-07-08"),as.Date("2019-07-15"),as.Date("2019-07-22"),
           as.Date("2019-07-29"),as.Date("2019-08-05"),as.Date("2019-08-12"),
           as.Date("2019-08-19"),as.Date("2019-08-28"),as.Date("2019-09-02"),
           as.Date("2019-09-11"),as.Date("2019-09-20"),as.Date("2019-09-27"),
           as.Date("2019-10-02"),as.Date("2019-10-11"),as.Date("2019-10-16"),
           as.Date("2019-10-23"),as.Date("2019-10-30"),as.Date("2019-11-07"))

# Sequence through the dates and the traps and execute the JAGS model
for(s in 1:length(dates)){
  
  # Select site and forecast date (#hashed out lines are for testing purposes)
  full_ebullition_model_alltrap_jags <- full_ebullition_model_alltrap %>% 
    filter(time <= dates[s])%>%
    #filter(time <= as.Date("2019-07-01"))%>%
    arrange(time)
  
  # fill in any missing covariate data This is done using impute TS --> but other models might be better
  full_ebullition_model_alltrap_jags$water_temp_dam <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags$water_temp_dam,option = "linear")
  full_ebullition_model_alltrap_jags$water_temp_dam_sd <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags$water_temp_dam_sd,option = "linear")
  full_ebullition_model_alltrap_jags$hobo_temp <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags$hobo_temp,option = "linear")
  full_ebullition_model_alltrap_jags$hobo_temp_sd <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags$hobo_temp_sd,option = "linear")
  full_ebullition_model_alltrap_jags$log_ebu_rate_sd <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags$log_ebu_rate_sd,option = "linear")
  full_ebullition_model_alltrap_jags$daily_temp <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags$daily_temp,option = "linear")
  full_ebullition_model_alltrap_jags$daily_temp_sd <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags$daily_temp_sd,option = "linear")
  full_ebullition_model_alltrap_jags$diff_bp <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags$diff_bp,option = "linear")
  full_ebullition_model_alltrap_jags$daily_wind <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags$daily_wind,option = "linear")
  full_ebullition_model_alltrap_jags$daily_wind_sd <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags$daily_wind_sd,option = "linear")
  
  # start the data frame used in JAGS on the first day ebullition was collected from the upstream traps
  y_nogaps <- full_ebullition_model_alltrap_jags$log_ebu_rate[!is.na(full_ebullition_model_alltrap_jags$log_ebu_rate)]
  sd_no_gaps <- full_ebullition_model_alltrap_jags$log_ebu_rate_sd[!is.na(full_ebullition_model_alltrap_jags$log_ebu_rate_sd)]
  
  # read the .nc file[h,s] from FLARE runs that correspond directly to dates ebullition from the traps were measured
  files <- list.files("./forecast_noaa_driver/")                                          
  start_forecast <- max(full_ebullition_model_alltrap_jags$time)         
  start <- paste(gsub("-", "", start_forecast))                  
  forecast_file <- grep(start, files, value=TRUE)
  noaa <- read_csv(paste0("./forecast_noaa_driver/",forecast_file))
  
  #prepare forecast temperatures from FLARE to be appended to the data frame recognized by jags
  forecast_jags <- noaa %>%
    select(forecast.date, ensembles, tmp2m, vgrd10m, ugrd10m, pressfc) %>%
    rename(time = forecast.date)%>% 
    mutate(time = as_date(time))%>%
    mutate(tmp2m = tmp2m - 273.15)%>%
    mutate(pressfc = pressfc * 0.001)%>%
    mutate(wind = sqrt(vgrd10m^2+ugrd10m^2))%>%
    group_by(time)%>%
    summarise(daily_bp = mean(pressfc),
              daily_bp_sd = sd(pressfc),
              daily_temp = mean(tmp2m),
              daily_temp_sd = sd(tmp2m),
              daily_wind = max(wind),
              daily_wind_sd = sd(wind))%>%
    mutate(diff_bp = daily_bp-lag(daily_bp))%>%
    select(time, daily_temp, daily_temp_sd, daily_wind, daily_wind_sd, diff_bp)%>%
    arrange(time)%>%
    mutate(hobo_temp = NA)%>%
    mutate(hobo_temp_sd = NA)%>%
    mutate(log_ebu_rate = NA)%>%
    mutate(log_ebu_rate_sd = NA)%>%
    mutate(water_temp_dam = NA)%>%
    mutate(water_temp_dam_sd = NA)%>%
    select(time, hobo_temp, hobo_temp_sd,daily_temp,daily_temp_sd,daily_wind,daily_wind_sd,diff_bp, water_temp_dam, water_temp_dam_sd, log_ebu_rate, log_ebu_rate_sd)%>%
    mutate(time = as.Date(time))%>%
    filter(time > start_forecast)
  
  
  # Prepare forecast so it can be saved as an .rds file and uploaded to github
  
  
  saveRDS(forecast_jags, paste0("./forecast_output/NOAA_forecast_alltraps",dates[s],".rds"))
  
  ### THIS IS THE ACTUAL JAGS MODELING PART ###
  # Append the temperature forecast from FLARE to the dataframe
  full_ebullition_model_alltrap_jags <- bind_rows(full_ebullition_model_alltrap_jags, forecast_jags)
  full_ebullition_model_alltrap_jags$hobo_temp_sd <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags$hobo_temp_sd,option = "linear")
  full_ebullition_model_alltrap_jags$log_ebu_rate_sd <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags$log_ebu_rate_sd,option = "linear")
  
 
  # Generate the data frame that is recognized by the jags model
  jags.data.ar = list(X = full_ebullition_model_alltrap_jags$log_ebu_rate, 
                      tau.obs = 1/(full_ebullition_model_alltrap_jags$log_ebu_rate_sd ^ 2),
                      N = nrow(full_ebullition_model_alltrap_jags), 
                      D_mean = full_ebullition_model_alltrap_jags$daily_temp,
                      tau.pre = 1/(full_ebullition_model_alltrap_jags$daily_temp_sd ^ 2),
                      B = full_ebullition_model_alltrap_jags$diff_bp)
  
  # Specify the parameters of interest from the model
  jags.params.ar = c("sd.pro", "mu2", "phi", "omega", "delta", "X", "Y", "D", "IC_forecast")
  
  #Initialize parameters 
  nchain = 3
  chain_seeds <- c(200,800,1400)
  init <- list()
  for(i in 1:nchain){
    init[[i]] <- list(sd.pro = sd(diff(y_nogaps)),
                      phi = rnorm(1,0.37, 0.09),
                      mu2 = rnorm(1, -5.11, 1.14),
                      omega = rnorm(1,0.30,0.08),
                      delta = rnorm(1,-0.53,0.25),
                      .RNG.name = "base::Wichmann-Hill", # so it is reproducible
                      .RNG.seed = chain_seeds[i])
  }
  
  #Initialize JAGS model
  j.model   <- jags.model(file = model.ar2,
                          data = jags.data.ar,
                          inits = init,
                          n.chains = 3)
  
  #Run JAGS model and sample from the posteriors
  eval_ebu  <- coda.samples(model = j.model,
                            variable.names = c("sd.pro", "mu2", "phi", "omega", "delta"),
                            n.iter = 10000, n.burnin = 1000)
  plot(eval_ebu)
  gelman_ebu <- gelman.diag(eval_ebu)
  gelman_ebu <- as.data.frame(bind_cols(gelman_ebu$mpsrf,as.Date(dates[s])))
  names(gelman_ebu) <- c("mpsrf","forecast_date")
  
  saveRDS(gelman_ebu, paste0("./forecast_output/noaa_ebu_model_gelman_diagnostics_",dates[s],".rds"))
  
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
  saveRDS(forecast_saved_ebu, paste0("./forecast_output/noaa_ebullition_predictive_interval_forecast_alltrap_",dates[s],".rds"))
  
  # Extract the parameter esimates from the ebullition model run for dates[s] and traps[h]
  ########################################################################################
  ebu_out_parms <- jags.out %>%
    spread_draws(sd.pro, mu2, phi, omega, delta, IC_forecast) %>%
    filter(.chain == 1) %>%
    rename(ensemble = .iteration) %>%
    mutate(forecast_date = start_forecast)%>%
    ungroup()%>%
    select(forecast_date, sd.pro, mu2, phi, omega,delta, IC_forecast)
  saveRDS(ebu_out_parms, paste0("./forecast_output/noaa_ebullition_parameters_",dates[s],".rds"))
  #######################################################################################
  
  
  #UNCERTAINTY PARTITIONING OF DATA
  ########################################################################################
  
  # initial condition
  
  hold_methane_pars <- T
  hold_methane_process <- T
  hold_FLARE_temp <- T
  hold_IC <- F
  
  # Here is the meat of the AR model that is using temperature and
  
  output <- matrix(ncol=420, nrow=16)
  
  for(m in 1:ncol(output)){
    
    parms <- sample_n(ebu_out_parms, 16)
    
    if(hold_methane_pars){
      mu2 <- sample(mean(ebu_out_parms$mu2, 16))
      phi <- sample(mean(ebu_out_parms$phi, 16))
      omega <- sample(mean(ebu_out_parms$omega, 16))
      delta <- sample(mean(ebu_out_parms$delta, 16))
    }else{
      mu2 <- as.vector(parms$mu2)
      phi <- as.vector(parms$phi)
      omega <- as.vector(parms$omega)
      delta <- as.vector(parms$delta)
    }
    
    if(hold_methane_process){
      process_error <- 0
    }else{
      process_error <- 1/parms$sd.pro^2
    }
    
    if(hold_FLARE_temp){
      temp <- forecast_jags$daily_temp
      pressure<- forecast_jags$diff_bp
    }else{
      temp <- rnorm(16, forecast_jags$daily_temp, forecast_jags$daily_temp_sd)
      pressure <- rnorm(16, forecast_jags$diff_bp, 0.1)
    }
    
    if(hold_IC){
      latent_ebu <- sample(mean(ebu_out_parms$IC_forecast, 16))
    }else{
      latent_ebu <- as.vector(parms$IC_forecast)
    }
    
    # This is the actual equation that is being run for the ebullition model
    
    output[,m] <- mu2 + phi * latent_ebu + omega * temp + delta*pressure + process_error
    ebu_forecast <- as.data.frame(output)
  }
  
  ebu_forecast <- cbind(forecast_saved_ebu$time[-1], ebu_forecast)
  names(ebu_forecast)[1] <- "time"
  ebu_forecast <- melt(ebu_forecast, id = "time")
  
  ebu_forecast_partition <- ebu_forecast %>%
    group_by(time) %>%
    summarise(mean = mean(value),
              max = max(value),
              min = min(value),
              upper = quantile(value, 0.90),
              lower = quantile(value, 0.10),
              var = var(value),
              sd = sd(value))%>%
    mutate(forecast_date = start_forecast)
  
  saveRDS(ebu_forecast_partition, paste0("./forecast_output/noaa_initial_condition_",dates[s],".rds"))
  
  #DRIVER DATA
  
  hold_methane_pars <- T
  hold_methane_process <- T
  hold_FLARE_temp <- F
  hold_IC <- T
  
  #Here is the meat of the AR model that is using temperature and
  
  output <- matrix(ncol=420, nrow=16)
  
  for(m in 1:ncol(output)){
    
    parms <- sample_n(ebu_out_parms, 16)
    
    if(hold_methane_pars){
      mu2 <- sample(mean(ebu_out_parms$mu2, 16))
      phi <- sample(mean(ebu_out_parms$phi, 16))
      omega <- sample(mean(ebu_out_parms$omega, 16))
      delta <- sample(mean(ebu_out_parms$delta, 16))
    }else{
      mu2 <- as.vector(parms$mu2)
      phi <- as.vector(parms$phi)
      omega <- as.vector(parms$omega)
      delta <- as.vector(parms$delta)
    }
    
    if(hold_methane_process){
      process_error <- 0
    }else{
      process_error <- 1/parms$sd.pro^2
    }
    
    if(hold_FLARE_temp){
      temp <- forecast_jags$daily_temp
      pressure<- forecast_jags$diff_bp
    }else{
      temp <- rnorm(16, forecast_jags$daily_temp, forecast_jags$daily_temp_sd)
      pressure <- rnorm(16, forecast_jags$diff_bp, 0.1)
    }
    
    if(hold_IC){
      latent_ebu <- sample(mean(ebu_out_parms$IC_forecast, 16))
    }else{
      latent_ebu <- as.vector(parms$IC_forecast)
    }
    
    # This is the actual equation that is being run for the ebullition model
    
    output[,m] <- mu2 + phi * latent_ebu + omega * temp + delta*pressure + process_error
    ebu_forecast <- as.data.frame(output)
  }
  
  ebu_forecast <- cbind(forecast_saved_ebu$time[-1], ebu_forecast)
  names(ebu_forecast)[1] <- "time"
  ebu_forecast <- melt(ebu_forecast, id = "time")
  
  ebu_forecast_partition <- ebu_forecast %>%
    group_by(time) %>%
    summarise(mean = mean(value),
              max = max(value),
              min = min(value),
              upper = quantile(value, 0.90),
              lower = quantile(value, 0.10),
              var = var(value),
              sd = sd(value))%>%
    mutate(forecast_date = start_forecast)
  
  saveRDS(ebu_forecast_partition, paste0("./forecast_output/noaa_driver_data_",dates[s],".rds"))
  
  #PROCESS
  
  hold_methane_pars <- T
  hold_methane_process <- F
  hold_FLARE_temp <- T
  hold_IC <- T
  
  #Here is the meat of the AR model that is using temperature and
  
  output <- matrix(ncol=420, nrow=16)
  
  for(m in 1:ncol(output)){
    
    parms <- sample_n(ebu_out_parms, 16)
    
    if(hold_methane_pars){
      mu2 <- sample(mean(ebu_out_parms$mu2, 16))
      phi <- sample(mean(ebu_out_parms$phi, 16))
      omega <- sample(mean(ebu_out_parms$omega, 16))
      delta <- sample(mean(ebu_out_parms$delta, 16))
    }else{
      mu2 <- as.vector(parms$mu2)
      phi <- as.vector(parms$phi)
      omega <- as.vector(parms$omega)
      delta <- as.vector(parms$delta)
    }
    
    if(hold_methane_process){
      process_error <- 0
    }else{
      process_error <- 1/parms$sd.pro^2
    }
    
    if(hold_FLARE_temp){
      temp <- forecast_jags$daily_temp
      pressure<- forecast_jags$diff_bp
    }else{
      temp <- rnorm(16, forecast_jags$daily_temp, forecast_jags$daily_temp_sd)
      pressure <- rnorm(16, forecast_jags$diff_bp, 0.1)
    }
    
    if(hold_IC){
      latent_ebu <- sample(mean(ebu_out_parms$IC_forecast, 16))
    }else{
      latent_ebu <- as.vector(parms$IC_forecast)
    }
    
    # This is the actual equation that is being run for the ebullition model
    
    output[,m] <- mu2 + phi * latent_ebu + omega * temp + delta*pressure + process_error
    ebu_forecast <- as.data.frame(output)
  }
  
  ebu_forecast <- cbind(forecast_saved_ebu$time[-1], ebu_forecast)
  names(ebu_forecast)[1] <- "time"
  ebu_forecast <- melt(ebu_forecast, id = "time")
  
  ebu_forecast_partition <- ebu_forecast %>%
    group_by(time) %>%
    summarise(mean = mean(value),
              max = max(value),
              min = min(value),
              upper = quantile(value, 0.90),
              lower = quantile(value, 0.10),
              var = var(value),
              sd = sd(value))%>%
    mutate(forecast_date = start_forecast)
  
  saveRDS(ebu_forecast_partition, paste0("./forecast_output/noaa_model_process_",dates[s],".rds"))
  
  #paramter
  
  hold_methane_pars <- F
  hold_methane_process <- T
  hold_FLARE_temp <- T
  hold_IC <- T
  
  #Here is the meat of the AR model that is using temperature and
  
  output <- matrix(ncol=420, nrow=16)
  
  for(m in 1:ncol(output)){
    
    parms <- sample_n(ebu_out_parms, 16)
    
    if(hold_methane_pars){
      mu2 <- sample(mean(ebu_out_parms$mu2, 16))
      phi <- sample(mean(ebu_out_parms$phi, 16))
      omega <- sample(mean(ebu_out_parms$omega, 16))
      delta <- sample(mean(ebu_out_parms$delta, 16))
    }else{
      mu2 <- as.vector(parms$mu2)
      phi <- as.vector(parms$phi)
      omega <- as.vector(parms$omega)
      delta <- as.vector(parms$delta)
    }
    
    if(hold_methane_process){
      process_error <- 0
    }else{
      process_error <- 1/parms$sd.pro^2
    }
    
    if(hold_FLARE_temp){
      temp <- forecast_jags$daily_temp
      pressure<- forecast_jags$diff_bp
    }else{
      temp <- rnorm(16, forecast_jags$daily_temp, forecast_jags$daily_temp_sd)
      pressure <- rnorm(16, forecast_jags$diff_bp, 0.1)
    }
    
    if(hold_IC){
      latent_ebu <- sample(mean(ebu_out_parms$IC_forecast, 16))
    }else{
      latent_ebu <- as.vector(parms$IC_forecast)
    }
    
    # This is the actual equation that is being run for the ebullition model
    
    output[,m] <- mu2 + phi * latent_ebu + omega * temp + delta*pressure + process_error
    ebu_forecast <- as.data.frame(output)
  }
  
  ebu_forecast <- cbind(forecast_saved_ebu$time[-1], ebu_forecast)
  names(ebu_forecast)[1] <- "time"
  ebu_forecast <- melt(ebu_forecast, id = "time")
  
  ebu_forecast_partition <- ebu_forecast %>%
    group_by(time) %>%
    summarise(mean = mean(value),
              max = max(value),
              min = min(value),
              upper = quantile(value, 0.90),
              lower = quantile(value, 0.10),
              var = var(value),
              sd = sd(value))%>%
    mutate(forecast_date = start_forecast)
  
  saveRDS(ebu_forecast_partition, paste0("./forecast_output/noaa_model_parameter_",dates[s],".rds"))
  
  ########################################################################################
}
