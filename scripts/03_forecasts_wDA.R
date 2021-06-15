#################################################################
# CH4cast                                                       #
# Ryan McClure                                                  #
# JAGS MODEL and FORECASTING SCRIPT                             #
#################################################################
set.seed(329)
#* TEMPERATURE SCALING MODEL ----
model.lm19 = ("lm19.txt")
jagsscript = cat("
model {  
   mu ~ dnorm(0,1e-6)  # intercet
   beta ~ dnorm(0,1e-6) # cat temp paramter
   sd.pro ~ dunif(0, 1000)
   tau.pro <-  pow(sd.pro, -2)
   
   for(i in 2:N) {
      predX[i] <- mu + C_mean[i]*beta 
      Y[i] ~ dnorm(predX[i],tau.pro) # model error
   }
}  
", 
                 file = model.lm19)


model.ar19 = ("ar19_model.txt")
jagsscript = cat("
model {  
   
   #priors===================================================
   
   mu2 ~ dnorm(0,1e-6)
   sd.pro ~ dunif(0, 1000)
   tau.pro <-  pow(sd.pro, -2)
   phi ~ dnorm(0,1e-6)
   omega ~ dnorm(0,1e-6)
   
   #Informative priors on initial conditions based on first observation
   X[1] ~ dnorm(x_init, tau.obs[1])
   
   #end priors===============================================
   
   for(i in 2:N) {
      
      #process model=============================================
      
      predX[i] <- mu2 + phi*X[i-1] + omega*D[i]
      X[i] ~ dnorm(predX[i],tau.pro)
      
      #end of process model======================================
      
      #data model================================================
      
      Y[i] ~ dnorm(X[i], tau.obs[i]) # Observation variation
            #end of data model=========================================
   }
  }", file = model.ar19)


# # Dates to forecast in 2019 --> Based off of dates starting from when the day ebullition was measured
dates <- c(as.Date("2019-06-17"),as.Date("2019-06-24"),as.Date("2019-07-01"),
           as.Date("2019-07-08"),as.Date("2019-07-15"),as.Date("2019-07-22"),
           as.Date("2019-07-29"),as.Date("2019-08-05"),as.Date("2019-08-12"),
           as.Date("2019-08-19"),as.Date("2019-08-28"),as.Date("2019-09-02"),
           as.Date("2019-09-11"),as.Date("2019-09-20"),as.Date("2019-09-27"),
           as.Date("2019-10-02"),as.Date("2019-10-11"),as.Date("2019-10-16"),
           as.Date("2019-10-23"),as.Date("2019-10-30"),as.Date("2019-11-07"))

 
t1 <- proc.time()

# Sequence through the dates and the traps and execute the JAGS model

for(s in 1:length(dates)){

    # Select site and forecast date (#hashed out lines are for testing purposes)
    full_ebullition_model_alltrap_jags <- full_ebullition_model_alltrap %>% 
      filter(time <= as.Date(dates[s]))%>%
      #filter(time <= as.Date("2019-07-01"))%>%
      #filter(time >= as.Date("2019-05-27"))%>%
      arrange(time)
    
    for (i in colnames(full_ebullition_model_alltrap_jags[,c(2:10,12)])) {
      full_ebullition_model_alltrap_jags[,i] <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags[,i],option = "linear")   
    }
    
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
  forecast_flare <- as.data.frame(temp[1:18,1:210,5:10])
  
  temp_forecast_1wk <- cbind(time, forecast_flare)%>%
    dplyr::filter(time >= dates[s])%>%
    dplyr::filter(time <= dates[s+1])%>%
    summarize_all(funs(mean))%>%
    select(-time)
  
  temp_forecast_2wk <- cbind(time, forecast_flare)%>%
    dplyr::filter(time > dates[s+1])%>%
    dplyr::filter(time <= dates[s+2])%>%
    summarize_all(funs(mean))%>%
    select(-time)
  
  temp_forecast <- dplyr::bind_rows(temp_forecast_1wk,temp_forecast_2wk)
    # Close out current NC file
    nc_close(nc)

    #* RUNJAGS FOR 2017 TEMP SCALE ----
    jags.data.lm = list(Y= full_ebullition_model_alltrap_jags$hobo_temp,
                        N = nrow(full_ebullition_model_alltrap_jags),
                        C_mean = full_ebullition_model_alltrap_jags$water_temp_dam)
    
    nchain = 3
    chain_seeds <- c(200,800,1400)
    init <- list()
    for(i in 1:nchain){
      init[[i]] <- list(sd.pro = temp_out_parms$mean_sd.pro,
                        mu = temp_out_parms$mean_mu,
                        beta = temp_out_parms$mean_beta,
                        .RNG.name = "base::Wichmann-Hill",
                        .RNG.seed = chain_seeds[i])
    }
    
    jags.params.lm.eval = c("sd.pro", "mu", "beta", "tau.pro")
    
    j.lm.model   <- jags.model(file = model.lm19,
                               data = jags.data.lm,
                               n.chains = 3)
    
    eval_temp  <- coda.samples(model = j.lm.model,
                               variable.names = jags.params.lm.eval,
                               n.iter = 10000, n.burnin = 1000)
    
    eval_temp  <- coda.samples(model = j.lm.model,
                               variable.names = jags.params.lm.eval,
                               n.iter = 30000, n.burnin = 20000)
    plot(eval_temp)
    
    # Extract the parameter esimates from the ebullition model run for dates[s] and traps[h]
    ########################################################################################
    ebu_temperature_parms <- eval_temp %>%
      spread_draws(sd.pro, mu, beta) %>%
      filter(.chain == 1) %>%
      rename(ensemble = .iteration) %>%
      mutate(forecast_date = start_forecast)%>%
      ungroup()%>%
      select(forecast_date, sd.pro, mu, beta)
    
    
    Nmc <- 1000         ## set number of Monte Carlo draws
    
    temp_forecast_function <- function(mu, beta, temp_forecast, Q, n = Nmc){
      est <- mu + (beta * temp_forecast) + Q
      return(est)
    }
    
    parms <- sample_n(ebu_temperature_parms, Nmc, replace=TRUE)
    
    # initial condition uncertainty
    forecast <- temp_forecast_function(temp_forecast = temp_forecast,
                                       mu = parms$mu,
                                       beta = parms$beta,
                                       Q = parms$sd.pro,
                                       n = Nmc)
    
    
      
    full_ebullition_model_alltrap_jags <- na.omit(full_ebullition_model_alltrap_jags)
    
    jags.data.ar = list(x_init = full_ebullition_model_alltrap_jags$log_ebu_rate[1],
                        Y = full_ebullition_model_alltrap_jags$log_ebu_rate, 
                        tau.obs = 1/((full_ebullition_model_alltrap_jags$log_ebu_rate_sd)/sqrt(4)) ^ 2,
                        N = nrow(full_ebullition_model_alltrap_jags), 
                        D = full_ebullition_model_alltrap_jags$hobo_temp)
    
    nchain = 3
    chain_seeds <- c(200,800,1400)
    init <- list()
    for(i in 1:nchain){
      init[[i]] <- list(sd.pro = ebu_out_parms$mean_sd.pro,
                        mu2 = ebu_out_parms$mean_mu2,
                        omega = ebu_out_parms$mean_omega, 
                        phi = ebu_out_parms$mean_phi,
                        .RNG.name = "base::Wichmann-Hill",
                        .RNG.seed = chain_seeds[i])
    }
    
    j.model   <- jags.model(file = model.ar19,
                            data = jags.data.ar,
                            inits = init,
                            n.chains = 3)
    
    jags.out  <- coda.samples(model = j.model,
                              variable.names = c("sd.pro", "mu2", "phi", "omega"),
                              n.iter = 30000, n.burnin = 10000)
    
    jags.out  <- coda.samples(model = j.model,
                              variable.names = c("sd.pro", "mu2", "phi", "omega"),
                              n.iter = 30000, n.burnin = 10000)

    plot(jags.out)
    gelman_ebu <- gelman.diag(jags.out)
    # gelman_ebu <- as.data.frame(bind_cols(gelman_ebu$mpsrf,as.Date(dates[s])))
    # names(gelman_ebu) <- c("mpsrf","forecast_date")
    # 
    # saveRDS(gelman_ebu, paste0("./forecast_output/ebu_model_gelman_diagnostics_",dates[s],".rds"))


    # Extract the parameter esimates from the ebullition model run for dates[s] and traps[h]
    ########################################################################################
    ebu_forecast_parms <- jags.out %>%
      spread_draws(sd.pro, mu2, phi, omega) %>%
      rename(ensemble = .iteration) %>%
      mutate(forecast_date = start_forecast)%>%
      ungroup()%>%
      select(forecast_date, sd.pro, mu2, phi, omega)
    saveRDS(ebu_out_parms, paste0("./forecast_output/ebullition_parameters_",dates[s],".rds"))
    #######################################################################################

    
    #UNCERTAINTY PARTITIONING OF DATA
    ########################################################################################

    Nmc <- 1000         ## set number of Monte Carlo draws
    future <- c(dates[s+1], dates[s+2])

    
    forecast_function <- function(IC, mu2, phi, omega, FLARE, Q, n = Nmc){
        est <- mu2 + phi * IC + omega * FLARE + Q
        return(est)
    }
    
    parms <- sample_n(ebu_forecast_parms, Nmc, replace=TRUE)
    
    # initial condition uncertainty
    forecast <- forecast_function(IC = rnorm(Nmc,full_ebullition_model_alltrap_jags$log_ebu_rate[s],full_ebullition_model_alltrap_jags$log_ebu_rate_sd[s]),  ## sample IC
                             FLARE = forecast,
                             mu2 = parms$mu2,
                             phi = parms$phi,
                             omega = parms$omega,
                             Q = parms$sd.pro,
                             n = Nmc)
    
    forecast <- cbind.data.frame(future, forecast)%>%
      rename(time = future)
    
    forecast <- melt(forecast, id = "time")
    
    forecast_ebu <- forecast %>%
      group_by(time) %>%
      summarise(mean = mean(value),
                max = max(value),
                min = min(value),
                upper = quantile(value, 0.95),
                lower = quantile(value, 0.05),
                var = var(value),
                sd = sd(value))%>%
      mutate(forecast_date = start_forecast)
    saveRDS(forecast_ebu, paste0("./forecast_output/ebu_forecast_",dates[s],".rds"))
    
    # # parameter uncertainty
    # param_uncertainty <- forecast_function(IC = mean(parms$IC_forecast),  ## sample IC
    #                                     FLARE = mean(FLARE$mean),
    #                                     mu2 = parms$mu2,
    #                                     phi = parms$phi,
    #                                     omega = parms$omega,
    #                                     Q = 1/sqrt(parms$sd.pro),
    #                                     n = Nmc)
    # 
    # param_uncertainty <- t(param_uncertainty)
    # param_uncertainty <- cbind.data.frame(future, param_uncertainty)%>%
    #   rename(time = future)
    # 
    # param_uncertainty <- melt(param_uncertainty, id = "time")
    # 
    # param_partition <- param_uncertainty %>%
    #   group_by(time) %>%
    #   summarise(mean = mean(value),
    #             max = max(value),
    #             min = min(value),
    #             upper = quantile(value, 0.95),
    #             lower = quantile(value, 0.05),
    #             var = var(value),
    #             sd = sd(value))%>%
    #   mutate(forecast_date = start_forecast)
    # saveRDS(param_partition, paste0("./forecast_output/model_parameter_",dates[s],".rds"))
    # 
    # 
    # # driver data uncertainty
    # driver_uncertainty <- forecast_function(IC = mean(parms$IC_forecast),  ## sample IC
    #                                        FLARE = rnorm(FLARE$mean,FLARE$sd),
    #                                        mu2 = mean(parms$mu2),
    #                                        phi = mean(parms$phi),
    #                                        omega = mean(parms$omega),
    #                                        Q = 0,
    #                                        n = Nmc)
    # 
    # driver_uncertainty <- t(driver_uncertainty)
    # driver_uncertainty <- cbind.data.frame(future, driver_uncertainty)%>%
    #   rename(time = future)
    # 
    # driver_uncertainty <- melt(driver_uncertainty, id = "time")
    # 
    # drive_partition <- driver_uncertainty %>%
    #   group_by(time) %>%
    #   summarise(mean = mean(value),
    #             max = max(value),
    #             min = min(value),
    #             upper = quantile(value, 0.95),
    #             lower = quantile(value, 0.05),
    #             var = var(value),
    #             sd = sd(value))%>%
    #   mutate(forecast_date = start_forecast)
    # saveRDS(drive_partition, paste0("./forecast_output/model_driver_",dates[s],".rds"))
    # 
    # # process uncertainty
    # process_uncertainty <- forecast_function(IC = mean(parms$IC_forecast),  ## sample IC
    #                                         FLARE = mean(FLARE$mean),
    #                                         mu2 = mean(parms$mu2),
    #                                         phi = mean(parms$phi),
    #                                         omega = mean(parms$omega),
    #                                         Q = 1/sqrt(parms$sd.pro),
    #                                         n = Nmc)
    # 
    # process_uncertainty <- t(process_uncertainty)
    # process_uncertainty <- cbind.data.frame(future, process_uncertainty)%>%
    #   rename(time = future)
    # 
    # process_uncertainty <- melt(process_uncertainty, id = "time")
    # 
    # process_partition <- process_uncertainty %>%
    #   group_by(time) %>%
    #   summarise(mean = mean(value),
    #             max = max(value),
    #             min = min(value),
    #             upper = quantile(value, 0.95),
    #             lower = quantile(value, 0.05),
    #             var = var(value),
    #             sd = sd(value))%>%
    #   mutate(forecast_date = start_forecast)
    # saveRDS(process_partition, paste0("./forecast_output/model_process_",dates[s],".rds"))
    
}

loop <- proc.time()-t1
print(loop[3]/60)
