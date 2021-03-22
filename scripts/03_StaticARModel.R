#################################################################
# CH4cast version 2                                             #
# Ryan McClure                                                  #
# JAGS MODEL and FORECASTING SCRIPT                             #
#################################################################

# The JAGS MODELs
model.lm = ("lm.txt")
jagsscript = cat("
model {  
   mu ~ dnorm(0,1);  # intercet
   beta ~ dnorm(0,1); # cat temp paramter
   sd.pro ~ dunif(0.00001, 10000);
   tau.pro <-  pow(sd.pro, -2)
   
   for(i in 2:N) {
      predX[i] <- mu + C[i]*beta; 
      X[i] ~ dnorm(predX[i],tau.pro); # Process variation
      Y[i] ~ dnorm(X[i], tau.obs[i]); # Observation variation
      C[i] ~ dnorm(C_mean[i],tau.pre[i]); # Covariate variation
   }
}  
", 
                 file = model.lm)

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
  
  # start the data frame used in JAGS on the first day ebullition was collected from the upstream traps
  y_nogaps <- full_ebullition_model_alltrap_jags$log_ebu_rate[!is.na(full_ebullition_model_alltrap_jags$log_ebu_rate)]
  sd_no_gaps <- full_ebullition_model_alltrap_jags$log_ebu_rate_sd[!is.na(full_ebullition_model_alltrap_jags$log_ebu_rate_sd)]
  
  # read the .nc file[h,s] from FLARE runs that correspond directly to dates ebullition from the traps were measured
  files <- list.files("./forecast_driver/")                                          
  start_forecast <- max(full_ebullition_model_alltrap_jags$time)         
  start <- paste(gsub("-", "_", start_forecast))                  
  start <-paste0(start,"_F_1")
  forecast_file <- grep(start, files, value=TRUE)
  nc <- nc_open(paste0("./forecast_driver/",forecast_file))
  t <- ncvar_get(nc,'time')
  time <- as.POSIXct(t, origin = '1970-01-01 00:00.00 UTC', tz = "EST")
  time <- strftime(time, format="%Y-%m-%d")
  temp <- ncvar_get(nc,'temp')
  
  #prepare forecast temperatures from FLARE to be appended to the data frame recognized by jags
  forecast_jags <- as.data.frame(cbind(time,temp[1:14,1:210,5], temp[1:14,1:210,6], temp[1:14,1:210,7],temp[1:14,1:210,8], temp[1:14,1:210,9], temp[1:14,1:210,10]))%>%
    filter(time>=start_forecast)%>%
    filter(time<=start_forecast+10)%>%
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
    select(time, hobo_temp, hobo_temp_sd, water_temp_dam, water_temp_dam_sd, log_ebu_rate, log_ebu_rate_sd)%>%
    mutate(time = as.Date(time))
  
  
  # Prepare forecast so it can be saved as an .rds file and uploaded to github
  forecast_save_temp <- as.data.frame(cbind(time,temp[1:14,1:210,5], temp[1:14,1:210,6], temp[1:14,1:210,7],temp[1:14,1:210,8], temp[1:14,1:210,9], temp[1:14,1:210,10]))%>%
    filter(time>=start_forecast)%>%
    filter(time<=start_forecast+10)%>%
    melt(., id = 'time')%>%
    group_by(time)%>%
    mutate(value = as.numeric(value))%>%
    summarise(mean = mean(value),
              max = max(value),
              min = min(value),
              upper = quantile(value, 0.95),
              lower = quantile(value, 0.05),
              var = var(value),
              sd = sd(value))
  
  saveRDS(forecast_save_temp, paste0("./forecast_output/FLARE_temp_forecast_alltraps",dates[s],".rds"))
  
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
  
  j.lm.model   <- jags.model(file = model.lm,
                             data = jags.data.lm,
                             n.chains = 3)
  
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

  forecast_time <- full_ebullition_model_alltrap_jags %>%
    filter(time >= start_forecast)%>%
    select(time)
  
  model_ic <- full_ebullition_model_alltrap_jags %>%
    filter(time == start_forecast)%>%
    select(log_ebu_rate, log_ebu_rate_sd)
  
  #RUN THE AR MODEL WITH PARAMETERS FROM 2017
  ########################################################################################
  
  # initial condition
  
  hold_methane_pars <- F
  hold_methane_process <- F
  hold_FLARE_temp <- F
  hold_IC <- F
  
  # Here is the meat of the AR model that is using temperature and
  
  output <- matrix(ncol=420, nrow=11)
  
  for(m in 1:ncol(output)){
    
    if(hold_methane_pars){
      mu2 <- sample(-6.58, 11, replace = T)
      phi <- sample(0.26, 11, replace = T)
      omega <- sample(0.38, 11, replace = T)
    }else{
      mu2 <- rnorm(11, -6.58, 1.27)
      phi <- rnorm(11, 0.26, 0.05)
      omega <- rnorm(11, 0.38, 0.05)
    }
    
    if(hold_methane_process){
      process_error <- 0
    }else{
      process_error <- 0
    }
    
    if(hold_FLARE_temp){
      FLARE_temp <- forecast_save_temp$mean
    }else{
      FLARE_temp <- rnorm(11, forecast_save_temp$mean, forecast_save_temp$sd)
    }
    
    if(hold_IC){
      latent_ebu <- model_ic$log_ebu_rate
    }else{
      latent_ebu <- rnorm(11, model_ic$log_ebu_rate, model_ic$log_ebu_rate_sd)
    }
    
    # This is the actual equation that is being run for the ebullition model
    
    output[,m] <- mu2 + phi * latent_ebu + omega * FLARE_temp + process_error
    ebu_forecast <- as.data.frame(output)
  }
  
  ebu_forecast <- cbind(forecast_time$time, ebu_forecast)
  names(ebu_forecast)[1] <- "time"
  ebu_forecast <- melt(ebu_forecast, id = "time")
  
  ebu_forecast_partition <- ebu_forecast %>%
    group_by(time) %>%
    summarise(mean = mean(value),
              upper_90 = quantile(value, 0.90),
              lower_90 = quantile(value, 0.10),
              upper_80 = quantile(value, 0.80),
              lower_80 = quantile(value, 0.20),
              upper_70 = quantile(value, 0.70),
              lower_70 = quantile(value, 0.30),
              upper_60 = quantile(value, 0.60),
              lower_60 = quantile(value, 0.40),
              var = var(value),
              sd = sd(value))%>%
    mutate(forecast_date = start_forecast)
  
  saveRDS(ebu_forecast_partition, paste0("./forecast_output/Static_AR_model_",dates[s],".rds"))
  
}
