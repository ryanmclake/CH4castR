#################################################################
# CH4cast version 2                                             #
# Ryan McClure                                                  #
# JAGS MODEL and FORECASTING SCRIPT woDA                        #
#################################################################


# # Dates to forecast in 2019 --> Based off of dates starting from when the day ebullition was measured

t1 <- proc.time()
# Sequence through the dates and the traps and execute the JAGS model
for(s in 1:length(dates)){
  
  # Select site and forecast date (#hashed out lines are for testing purposes)
  full_ebullition_model_alltrap_jags <- full_ebullition_model_alltrap %>% 
    filter(time <= as.Date(dates[s]))%>%
    filter(time >= as.Date("2019-05-27"))%>%
    mutate(hobo_temp = ifelse(time<=as.Date("2019-06-17"),rnorm(length(.),mean = mean(full_ebullition_model_17[21:42,2]),sd = 1),hobo_temp),
           hobo_temp_sd = ifelse(time<=as.Date("2019-06-17"),rnorm(length(.),mean = mean(full_ebullition_model_17[21:42,3]),sd = 0.5),hobo_temp_sd))%>%
    arrange(time)
  
  # fill in any missing covariate data This is done using impute TS --> but other models might be better
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
  
  jags.params.lm.eval = c("sd.pro", "mu", "beta")
  
  
  j.lm.model   <- jags.model(file = model.lm19,
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
                      phi = rnorm(1,0.26, 0.05),
                      mu2 = rnorm(1, -6.58, 1.27),
                      omega = rnorm(1,0.38,0.05),
                      .RNG.name = "base::Wichmann-Hill", # so it is reproducible
                      .RNG.seed = chain_seeds[i])
  }
  
  #Initialize JAGS model
  j.model   <- jags.model(file = model.ar19,
                          data = jags.data.ar,
                          inits = init,
                          n.chains = 3)
  
  
  #Run JAGS model and sample from the posteriors
  jags.out   <- coda.samples(model = j.model,
                             variable.names = jags.params.ar, 
                             n.iter = 10000, n.burnin = 1000)

  
  
  # Extract the parameter esimates from the ebullition model run for dates[s] and traps[h]
  ########################################################################################
  ebu_out_parms <- jags.out %>%
    spread_draws(sd.pro, mu2, phi, omega, IC_forecast) %>%
    filter(.chain == 1) %>%
    rename(ensemble = .iteration) %>%
    mutate(forecast_date = start_forecast)%>%
    ungroup()%>%
    select(forecast_date, sd.pro, mu2, phi, omega, IC_forecast)
  #######################################################################################
  
  
  #UNCERTAINTY PARTITIONING OF DATA
  ########################################################################################
  
  
  Nmc <- 1000         ## set number of Monte Carlo draws
  future <- seq(from = start_forecast, to = start_forecast+lubridate::days(16), by = "days")
  ylim <- c(-10, 10)  ## set Y range on plot
  N.cols <- c("black", "red", "green", "blue", "orange") ## set colors
  trans <- 0.8       ## set transparancy
  
  
  forecast_function <- function(IC, mu2, phi, omega, FLARE, Q = 0, n = Nmc){
    N <- matrix(NA, n, 17)  ## storage
    Npre <- IC
    for(t in 1:ncol(N)){
      est <- mu2 + phi * Npre + omega * FLARE + Q
      N[,t] <- rnorm(n, est, Q)                                    ## predict next step
      Npre <- N[, t]                                              ## update IC
    }
    return(N)
  }
  
  parms <- sample_n(ebu_out_parms, Nmc, replace=TRUE)
  FLARE = scale_out_forecast[,c(1,2,12)]%>%filter(time>=start_forecast)
  drow <- sample.int(nrow(forecast_temp_partition), Nmc, replace = TRUE)

  
  # process uncertainty
  no_da_uncertainty <- forecast_function(IC = parms$IC_forecast,  ## sample IC
                                           FLARE = FLARE$mean,
                                           mu2 = -6.58,
                                           phi = 0.26,
                                           omega = 0.38,
                                           Q = 1/sqrt(parms$sd.pro),
                                           n = Nmc)
  
  no_da_uncertainty <- t(no_da_uncertainty)
  no_da_uncertainty <- cbind.data.frame(future, no_da_uncertainty)%>%
    rename(time = future)
  
  no_da_uncertainty <- melt(no_da_uncertainty, id = "time")
  
  no_da_partition <- no_da_uncertainty %>%
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
              sd = sd(value),.groups = "drop")%>%
    mutate(forecast_date = start_forecast)
  saveRDS(no_da_partition, paste0("./forecast_output/no_da_forecast_",dates[s],".rds"))
}

dopar_loop <- proc.time()-t1
print(dopar_loop[3]/60)