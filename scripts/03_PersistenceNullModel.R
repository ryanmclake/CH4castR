#################################################################
# CH4cast version 2                                             #
# Ryan McClure                                                  #
# PERSISTENCE NULL MODEL --> BASED FROM ECO4CAST CHALLENGE      #
#################################################################

set.seed(329)

# number of days to forecast --> based on forecast horizon from the .nc files from FLARE
f_days = 10

# Dates to forecast
dates <- c(as.Date("2019-06-17"),as.Date("2019-06-24"),as.Date("2019-07-01"),
           as.Date("2019-07-08"),as.Date("2019-07-15"),as.Date("2019-07-22"),
           as.Date("2019-07-29"),as.Date("2019-08-05"),as.Date("2019-08-12"),
           as.Date("2019-08-19"),as.Date("2019-08-28"),as.Date("2019-09-02"),
           as.Date("2019-09-11"),as.Date("2019-09-20"),as.Date("2019-09-27"),
           as.Date("2019-10-02"),as.Date("2019-10-11"),as.Date("2019-10-16"),
           as.Date("2019-10-23"),as.Date("2019-10-30"),as.Date("2019-11-07"))

#'Generic random walk state-space model is JAGS format.  We use this model for 
#'both the oxygen and temperature null forecasts
RandomWalk = "
model{
  # Priors
  x[1] ~ dnorm(x_ic,tau_init)
  tau_add ~ dgamma(0.1,0.1)
  tau_init ~ dgamma(0.1,0.1)
  
  # Process Model
  for(t in 2:n){
    x[t]~dnorm(x[t-1],tau_add)
    x_obs[t] ~ dnorm(x[t],tau_obs[t])
  }
  # Data Model
  for(i in 1:nobs){
    y[i] ~ dnorm(x[y_wgaps_index[i]], tau_obs[y_wgaps_index[i]])
  }
}
"

#'Create variable for combined forecasts across sites
for(s in 1:length(dates)){
  
    targets <- full_ebullition_model_alltrap %>% 
      #filter(time <= "2019-07-01")
      filter(time <= dates[s])
    
    targets$water_temp_dam <- imputeTS::na_interpolation(targets$water_temp_dam,option = "linear")
    targets$water_temp_dam_sd <- imputeTS::na_interpolation(targets$water_temp_dam_sd,option = "linear")
    targets[1,7] <- 0.5
    targets$log_ebu_rate_sd <- imputeTS::na_interpolation(targets$log_ebu_rate_sd,option = "linear")
    
    
  # Select site
  site_data_var <- targets
  # Find the last day in the observed data and add one day for the start of the 
  # forecast
  start_forecast <- max(site_data_var$time)
  
  # This is key here - I added the forecast horizon on the end of the data for the forecast period
  full_time <- tibble(time = seq(min(site_data_var$time), max(site_data_var$time) + days(f_days), by = "1 day"))
  
  # Join the full time with the site_data_var so there aren't gaps in the time column
  site_data_var <- left_join(full_time, site_data_var)
  
  #observed ebullition: Full time series with gaps
  y_wgaps <- site_data_var$log_ebu_rate
  sd_wgaps <- imputeTS::na_interpolation(site_data_var$log_ebu_rate_sd,option = "linear")
  time <- c(site_data_var$time)
  #observed oxygen: time series without gaps
  y_nogaps <- y_wgaps[!is.na(y_wgaps)]
  #Index: time series with gaps
  y_wgaps_index <- 1:length(y_wgaps)
  #Index: the index of the non-NA values in time series with gaps
  y_wgaps_index <- y_wgaps_index[!is.na(y_wgaps)]
  
  #Generate starting initial conditions for latent states
  init_x <- approx(x = time[!is.na(y_wgaps)], y = y_nogaps, xout = time, rule = 2)$y
  
  #Create a list of the data for use in JAGS.  Include vector lengths (nobs, n)
  data <- list(y = y_nogaps,
               y_wgaps_index = y_wgaps_index,
               nobs = length(y_wgaps_index),
               tau_obs = 1/(sd_wgaps ^ 2),
               n = length(y_wgaps),
               x_ic = 0.0)
  #Initialize parameters 
  nchain = 3
  chain_seeds <- c(200,800,1400)
  init <- list()
  for(i in 1:nchain){
    init[[i]] <- list(tau_add = 1/var(diff(y_nogaps)),
                      tau_init = mean( 1/var(diff(y_nogaps)), na.rm = TRUE),
                      .RNG.name = "base::Wichmann-Hill",
                      .RNG.seed = chain_seeds[i],
                      x = init_x)
  }
  
  #Initialize JAGS model
  j.model   <- jags.model (file = textConnection(RandomWalk),
                           data = data,
                           inits = init,
                           n.chains = 3)
  
  
  #Run JAGS model as the burn-in
  jags.out   <- coda.samples(model = j.model,variable.names = c("tau_add","tau_init"), n.iter = 10000)
  
  #Run JAGS model again and sample from the posteriors
  m   <- coda.samples(model = j.model,
                      variable.names = c("x","tau_add","tau_init", "x_obs"),
                      n.iter = 10000,
                      thin = 5)
  
  #Use TidyBayes package to clean up the JAGS output
  model_output <- m %>%
    spread_draws(x_obs[day]) %>%
    filter(.chain == 1) %>%
    rename(ensemble = .iteration) %>%
    mutate(time = full_time$time[day]) %>%
    ungroup() %>%
    select(time, x_obs, ensemble)%>%
    group_by(time) %>% 
    summarise(mean = mean(x_obs),
              upper_90 = quantile(x_obs, 0.90),
              lower_90 = quantile(x_obs, 0.10),
              upper_80 = quantile(x_obs, 0.80),
              lower_80 = quantile(x_obs, 0.20),
              upper_70 = quantile(x_obs, 0.70),
              lower_70 = quantile(x_obs, 0.30),
              upper_60 = quantile(x_obs, 0.60),
              lower_60 = quantile(x_obs, 0.40),
              var = var(x_obs),
              sd = sd(x_obs))%>%
    filter(time>=start_forecast)%>%
    mutate(forecast_date = start_forecast)
  saveRDS(model_output, paste0("./forecast_output/ebullition_null_persistence_forecast_alltrap_",dates[s],".rds"))
}

