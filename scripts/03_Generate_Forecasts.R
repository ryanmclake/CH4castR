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
   
   pars ~ dmnorm(prior_mean,prior_inv_cov)
   sd.pro ~ dlnorm(proc_mean, proc_prec)
   
   for(i in 2:N){
      tau.pro[i] <- 1/((sd.pro*ndays[i])*(sd.pro*ndays[i]))
      predX[i] <- pars[1] + C_mean[i]*pars[2] 
      Y[i] ~ dnorm(predX[i],tau.pro[i])
   }
}  
", file = model.lm19)


model.ar19 = ("ar19_model.txt")
jagsscript = cat("
model {  
   
   #priors===================================================
   pars ~ dmnorm(prior_mean,prior_inv_cov)
   sd.pro ~ dlnorm(proc_mean, proc_prec)
   
   #Informative priors on initial conditions based on first observation
   X[1] ~ dnorm(x_init, tau.obs[1])
   
   #end priors===============================================
   
   for(i in 2:N) {
      
      #process model=============================================
      
      tau.pro[i] <- 1/((sd.pro*ndays[i])*(sd.pro*ndays[i]))
      predX[i] <- pars[1] + pars[2]*X[i-1] + pars[3]*D[i]
      X[i] ~ dnorm(predX[i],tau.pro[i])
      
      #end of process model======================================
      
      #data model================================================
      
      Y[i] ~ dnorm(X[i], tau.obs[i]) # Observation variation
            #end of data model=========================================
   }
  }", file = model.ar19)


null = ("null_model.txt")
jagsscript = cat("
model {  
   
   #priors===================================================
   
   sd.pro ~ dlnorm(proc_mean, proc_prec)
   
   #Informative priors on initial conditions based on first observation
   X[1] ~ dnorm(x_init, tau.obs[1])
   
   #end priors===============================================
   
   for(i in 2:N) {
      
      #process model=============================================
      
      tau.pro[i] <- 1/((sd.pro*ndays[i])*(sd.pro*ndays[i]))
      predX[i] <- X[i-1]
      X[i] ~ dnorm(predX[i],tau.pro[i])
      
      
      #end of process model======================================
      
      #data model================================================
      
      Y[i] ~ dnorm(X[i], tau.obs[i]) # Observation variation
            #end of data model=========================================
   }
  }", file = null)


# # Dates to forecast in 2019 --> Based off of dates starting from when the day ebullition was measured
dates <- c(as.Date("2019-06-17"),as.Date("2019-06-24"),as.Date("2019-07-01"),
           as.Date("2019-07-08"),as.Date("2019-07-15"),as.Date("2019-07-22"),
           as.Date("2019-07-29"),as.Date("2019-08-05"),as.Date("2019-08-12"),
           as.Date("2019-08-19"),as.Date("2019-08-28"),as.Date("2019-09-02"),
           as.Date("2019-09-11"),as.Date("2019-09-20"),as.Date("2019-09-27"),
           as.Date("2019-10-02"),as.Date("2019-10-11"),as.Date("2019-10-16"),
           as.Date("2019-10-23"),as.Date("2019-10-30"),as.Date("2019-11-07"),
           as.Date("2019-11-14"),as.Date("2019-11-21"))


t1 <- proc.time()

# Sequence through the dates and the traps and execute the JAGS model

for(s in 1:length(dates)){
  
  if (dates[s] == as.Date("2019-11-14")) break
  
  full_ebullition_model_alltrap_jags <- full_ebullition_model_alltrap%>%
    select(time, ebu_rate, ebu_rate_se, hobo_temp, water_temp_dam)%>%
    filter(time <= dates[s])%>%
    filter(time >= "2019-05-27")
  
  for (i in colnames(full_ebullition_model_alltrap_jags[,c(3:5)])) {
    full_ebullition_model_alltrap_jags[,i] <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags[,i],option = "linear")   
  }
  
  full_ebullition_model_alltrap_jags <- na.omit(full_ebullition_model_alltrap_jags)%>%
    mutate(ndays = difftime(time,lag(time)))
  
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
  forecast_flare3 <- cbind(time,as.data.frame(temp[1:18,1:210,7]))
  forecast_flare6 <- cbind(time,as.data.frame(temp[1:18,1:210,10]))
  
  forecast_flare <- bind_rows(forecast_flare3,forecast_flare6) %>%
    group_by(time) %>%
    summarise_all(funs(mean))%>%
    select(-time)
  
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
                      C_mean = full_ebullition_model_alltrap_jags$water_temp_dam,
                      ndays = full_ebullition_model_alltrap_jags$ndays,
                      prior_mean = colMeans(temp_out_parms[,5:6]),
                      prior_inv_cov = solve(cov(temp_out_parms[,5:6])),
                      proc_mean = log((mean(temp_out_parms$sd.pro)^2)/sqrt(mean(temp_out_parms$sd.pro)^2 + var(temp_out_parms$sd.pro))),
                      proc_prec = 1 / log(1 + (var(temp_out_parms$sd.pro)/(mean(temp_out_parms$sd.pro)^2))))
  
  nchain = 3
  chain_seeds <- c(200,800,1400)
  init <- list()
  for(i in 1:nchain){
    init[[i]] <- list(sd.pro = runif(1, 0.5,1),
                      pars = c(runif(1,1,4),runif(1, 0.6, 0.8)),
                      .RNG.name = "base::Wichmann-Hill",
                      .RNG.seed = chain_seeds[i])
  }
  
  jags.params.lm.eval = c("sd.pro", "pars[1]", "pars[2]")
  
  j.lm.model   <- jags.model(file = model.lm19,
                             data = jags.data.lm,
                             n.chains = 3)
  
  eval_temp  <- coda.samples(model = j.lm.model,
                             variable.names = jags.params.lm.eval,
                             n.iter = 10000, n.burnin = 1000)
  
  
  #plot(eval_temp)
  
  # Extract the parameter esimates from the ebullition model run for dates[s] and traps[h]
  ########################################################################################
  temperature_scale_parms <- eval_temp %>%
    spread_draws(sd.pro, `pars[1]`, `pars[2]`) %>%
    filter(.chain == 1) %>%
    rename(ensemble = .iteration) %>%
    mutate(forecast_date = start_forecast)%>%
    ungroup()%>%
    select(forecast_date, sd.pro, `pars[1]`, `pars[2]`)
  
  
  temp_forecast_function <- function(mu, beta, temp_forecast, Q, ndays){
    est <- mu + (beta * temp_forecast) + rnorm(nrow(full_ebullition_model_alltrap_jags),0, sd = Q*ndays)
    return(est)
  }
  
  parms <- sample_n(temperature_scale_parms, 210, replace=TRUE)
  
  ndays = full_ebullition_model_alltrap_jags %>% filter(time == dates[s])%>%
    select(ndays)
  
  # initial condition uncertainty
  forecast_1wk <- temp_forecast_function(temp_forecast = temp_forecast[1,],
                                         mu = parms$`pars[1]`,
                                         beta = parms$`pars[2]`,
                                         Q = parms$sd.pro,
                                         ndays = as.numeric(ndays))
  
  forecast_2wk <- temp_forecast_function(temp_forecast = temp_forecast[2,],
                                         mu = parms$`pars[1]`,
                                         beta = parms$`pars[2]`,
                                         Q = parms$sd.pro,
                                         ndays = as.numeric(ndays))
  
  jags.data.ar = list(x_init = full_ebullition_model_alltrap_jags$ebu_rate[1],
                      Y = full_ebullition_model_alltrap_jags$ebu_rate, 
                      tau.obs = 1/((full_ebullition_model_alltrap_jags$ebu_rate_se)) ^ 2,
                      N = nrow(full_ebullition_model_alltrap_jags), 
                      D = full_ebullition_model_alltrap_jags$hobo_temp,
                      ndays = full_ebullition_model_alltrap_jags$ndays,
                      prior_mean = colMeans(ebu_out_parms[,5:7]),
                      prior_inv_cov = solve(cov(ebu_out_parms[,5:7])),
                      proc_prec = 1 / log(1 + (var(ebu_out_parms$sd.pro)/(mean(ebu_out_parms$sd.pro)^2))),
                      proc_mean = log((mean(ebu_out_parms$sd.pro)^2)/sqrt(mean(ebu_out_parms$sd.pro)^2 + var(ebu_out_parms$sd.pro))))
  
  nchain = 3
  chain_seeds <- c(200,800,1400)
  init <- list()
  for(i in 1:nchain){
    init[[i]] <- list(sd.pro = runif(1, 0.01, 2),
                      mu2 = runif(1, -20,0),
                      omega = runif(1, 0.3,1), 
                      phi = runif(1, -0.2, 0.3),
                      .RNG.name = "base::Wichmann-Hill",
                      .RNG.seed = chain_seeds[i])
  }
  
  j.model   <- jags.model(file = model.ar19,
                          data = jags.data.ar,
                          inits = init,
                          n.chains = 3)
  
  jags.out  <- coda.samples(model = j.model,
                            variable.names = c("sd.pro", "pars[1]", "pars[2]", "pars[3]"),
                            n.iter = 200000, n.burnin = 20000, thin = 200)
  
  gelman_ebu <- gelman.diag(jags.out)
  gelman_ebu <- as.data.frame(bind_cols(gelman_ebu$mpsrf,as.Date(dates[s])))
  names(gelman_ebu) <- c("mpsrf","forecast_date")
  
  saveRDS(gelman_ebu, paste0("./forecast_output/ebu_model_gelman_diagnostics_",dates[s],".rds"))
  
  
  # Extract the parameter esimates from the ebullition model run for dates[s] and traps[h]
  ########################################################################################
  ebu_forecast_parms <- jags.out %>%
    spread_draws(sd.pro, `pars[1]`, `pars[2]`, `pars[3]`) %>%
    rename(ensemble = .iteration) %>%
    mutate(forecast_date = start_forecast)%>%
    ungroup()%>%
    select(forecast_date, sd.pro, `pars[1]`, `pars[2]`, `pars[3]`)
  saveRDS(ebu_forecast_parms, paste0("./forecast_output/ebullition_parameters_",dates[s],".rds"))
  #######################################################################################
  
  
  #UNCERTAINTY PARTITIONING OF DATA
  ########################################################################################
  
  ## set number of Monte Carlo draws
  time <- c(dates[s+1], dates[s+2])
  
  
  forecast_function <- function(IC, mu2, phi, omega, FLARE, Q, ndays){
    est <- mu2 + (phi * IC) + (omega * FLARE) + rnorm(210, 0, sd = Q*ndays)
    return(est)
  }
  
  
  parms <- sample_n(ebu_forecast_parms, 210, replace=TRUE)
  
  IC <- full_ebullition_model_alltrap_jags %>% filter(time == dates[s]) %>%
    select(ebu_rate, ebu_rate_se)
  
  
  observed_date <- full_ebullition_model_alltrap_jags %>% filter(time == dates[s]) %>%
    select(time)
  
  # initial condition uncertainty
  ebu_forecast_1wk <- forecast_function(IC = rnorm(210,IC[,1], IC[,2]), ## sample IC
                                        FLARE = forecast_1wk,
                                        mu2 = parms$`pars[1]`,
                                        phi = parms$`pars[2]`,
                                        omega = parms$`pars[3]`,
                                        Q = parms$sd.pro,
                                        ndays = as.numeric(ndays))
  
  ebu_forecast_2wk <- forecast_function(IC = ebu_forecast_1wk,  ## sample IC
                                        FLARE = forecast_2wk,
                                        mu2 = parms$`pars[1]`,
                                        phi = parms$`pars[2]`,
                                        omega = parms$`pars[3]`,
                                        Q = parms$sd.pro,
                                        ndays = as.numeric(ndays))
  
  ebu_forecast <- rbind(ebu_forecast_1wk, ebu_forecast_2wk)%>% unname(.)%>%
    cbind(time,.)
  
  
  observed <- cbind(observed_date, t(rnorm(210,IC[,1], IC[,2])))
  
  ebu_forecast <- bind_rows(observed, ebu_forecast)
  
  ebu_forecast_raw <- melt(ebu_forecast, id = "time")
  
  final_forecast_sumarized <- ebu_forecast_raw %>%
    group_by(time) %>%
    summarise(mean = mean(value),
              max = max(value),
              min = min(value),
              upper_95 = quantile(value, 0.95, na.rm = T),
              lower_95 = quantile(value, 0.05, na.rm = T),
              var = var(value))%>%
    mutate(forecast_date = start_forecast)
  
  saveRDS(final_forecast_sumarized, paste0("./forecast_output/ebu_forecast_wDA_",dates[s],".rds"))
  
  ### BELOW partitions uncertainty ###
  # parameter uncertainty
  parm_ebu_forecast_1wk <- forecast_function(IC = mean(rnorm(210,IC[,1], IC[,2])), ## sample IC
                                             FLARE = rowMeans(forecast_1wk),
                                             mu2 = parms$`pars[1]`,
                                             phi = parms$`pars[2]`,
                                             omega = parms$`pars[3]`,
                                             Q = 0,
                                             ndays = as.numeric(ndays))
  
  parm_ebu_forecast_2wk <- forecast_function(IC = mean(parm_ebu_forecast_1wk), ## sample IC
                                             FLARE = rowMeans(forecast_2wk),
                                             mu2 = parms$`pars[1]`,
                                             phi = parms$`pars[2]`,
                                             omega = parms$`pars[3]`,
                                             Q = 0,
                                             ndays = as.numeric(ndays))
  
  parm_ebu_forecast <- rbind(parm_ebu_forecast_1wk, parm_ebu_forecast_2wk)%>% unname(.)%>%
    cbind.data.frame(time,.)
  
  
  observed <- cbind(observed_date, t(rnorm(210,IC[,1], IC[,2])))
  
  parm_ebu_forecast <- bind_rows(observed, parm_ebu_forecast)
  
  parm_ebu_forecast_raw <- melt(parm_ebu_forecast, id = "time")
  
  parm_final_forecast_sumarized <- parm_ebu_forecast_raw %>%
    group_by(time) %>%
    summarise(mean = mean(value),
              max = max(value),
              min = min(value),
              upper_95 = quantile(value, 0.95, na.rm = T),
              lower_95 = quantile(value, 0.05, na.rm = T),
              var = var(value))%>%
    mutate(forecast_date = start_forecast)
  
  saveRDS(parm_final_forecast_sumarized, paste0("./forecast_output/model_parameter",dates[s],".rds"))
  
  
  # driver data uncertainty
  driv_ebu_forecast_1wk <- forecast_function(IC = mean(rnorm(210,IC[,1], IC[,2])), ## sample IC
                                             FLARE = forecast_1wk,
                                             mu2 = mean(parms$`pars[1]`),
                                             phi = mean(parms$`pars[2]`),
                                             omega = mean(parms$`pars[3]`),
                                             Q = 0,
                                             ndays = as.numeric(ndays))
  
  driv_ebu_forecast_2wk <- forecast_function(IC = rowMeans(driv_ebu_forecast_1wk), ## sample IC
                                             FLARE = forecast_2wk,
                                             mu2 = mean(parms$`pars[1]`),
                                             phi = mean(parms$`pars[2]`),
                                             omega = mean(parms$`pars[3]`),
                                             Q = 0,
                                             ndays = as.numeric(ndays))
  
  driv_ebu_forecast <- rbind(driv_ebu_forecast_1wk, driv_ebu_forecast_2wk)%>% unname(.)%>%
    cbind.data.frame(time,.)
  
  
  observed <- cbind(observed_date, t(rnorm(210,IC[,1], IC[,2])))
  
  driv_ebu_forecast <- bind_rows(observed, driv_ebu_forecast)
  
  driv_ebu_forecast_raw <- melt(driv_ebu_forecast, id = "time")
  
  driv_final_forecast_sumarized <- driv_ebu_forecast_raw %>%
    group_by(time) %>%
    summarise(mean = mean(value),
              max = max(value),
              min = min(value),
              upper_95 = quantile(value, 0.95, na.rm = T),
              lower_95 = quantile(value, 0.05, na.rm = T),
              var = var(value))%>%
    mutate(forecast_date = start_forecast)
  saveRDS(driv_final_forecast_sumarized, paste0("./forecast_output/model_driver_",dates[s],".rds"))
  
  # process uncertainty
  proc_ebu_forecast_1wk <- forecast_function(IC = mean(rnorm(210,IC[,1], IC[,2])), ## sample IC
                                             FLARE = as.numeric(rowMeans(forecast_1wk)),
                                             mu2 = mean(parms$`pars[1]`),
                                             phi = mean(parms$`pars[2]`),
                                             omega = mean(parms$`pars[3]`),
                                             Q = parms$sd.pro,
                                             ndays = as.numeric(ndays))
  
  proc_ebu_forecast_2wk <- forecast_function(IC = mean(proc_ebu_forecast_1wk), ## sample IC
                                             FLARE = as.numeric(rowMeans(forecast_2wk)),
                                             mu2 = mean(parms$`pars[1]`),
                                             phi = mean(parms$`pars[2]`),
                                             omega = mean(parms$`pars[3]`),
                                             Q = parms$sd.pro,
                                             ndays = as.numeric(ndays))
  
  proc_ebu_forecast <- rbind(proc_ebu_forecast_1wk, proc_ebu_forecast_2wk)%>% unname(.)%>%
    cbind.data.frame(time,.)
  
  
  observed <- cbind(observed_date, t(rnorm(210,IC[,1], IC[,2])))
  
  proc_ebu_forecast <- bind_rows(observed, proc_ebu_forecast)
  
  proc_ebu_forecast_raw <- melt(proc_ebu_forecast, id = "time")
  
  proc_ebu_forecast_sumarized <- proc_ebu_forecast_raw %>%
    group_by(time) %>%
    summarise(mean = mean(value),
              max = max(value),
              min = min(value),
              upper_95 = quantile(value, 0.95, na.rm = T),
              lower_95 = quantile(value, 0.05, na.rm = T),
              var = var(value))%>%
    mutate(forecast_date = start_forecast)
  saveRDS(proc_ebu_forecast_sumarized, paste0("./forecast_output/model_process_",dates[s],".rds"))
  
  
  # initial conditions
  inic_ebu_forecast_1wk <- forecast_function(IC = rnorm(210,IC[,1], IC[,2]), ## sample IC
                                             FLARE = rowMeans(forecast_1wk),
                                             mu2 = mean(parms$`pars[1]`),
                                             phi = mean(parms$`pars[2]`),
                                             omega = mean(parms$`pars[3]`),
                                             Q = 0,
                                             ndays = as.numeric(ndays))
  
  inic_ebu_forecast_2wk <- forecast_function(IC = inic_ebu_forecast_1wk, ## sample IC
                                             FLARE = rowMeans(forecast_2wk),
                                             mu2 = mean(parms$`pars[1]`),
                                             phi = mean(parms$`pars[2]`),
                                             omega = mean(parms$`pars[3]`),
                                             Q = 0,
                                             ndays = as.numeric(ndays))
  
  inic_ebu_forecast <- rbind(inic_ebu_forecast_1wk, inic_ebu_forecast_2wk)%>% unname(.)%>%
    cbind.data.frame(time,.)
  
  
  observed <- cbind(observed_date, t(rnorm(210,IC[,1], IC[,2])))
  
  inic_ebu_forecast <- bind_rows(observed, inic_ebu_forecast)
  
  inic_ebu_forecast_raw <- melt(inic_ebu_forecast, id = "time")
  
  inic_ebu_forecast_sumarized <- inic_ebu_forecast_raw %>%
    group_by(time) %>%
    summarise(mean = mean(value),
              max = max(value),
              min = min(value),
              upper_95 = quantile(value, 0.95, na.rm = T),
              lower_95 = quantile(value, 0.05, na.rm = T),
              var = var(value))%>%
    mutate(forecast_date = start_forecast)
  saveRDS(inic_ebu_forecast_sumarized, paste0("./forecast_output/model_initial_",dates[s],".rds"))
  
  
  ### Generate the forecasts without data assimilation
  
  parm_17 <- sample_n(ebu_out_parms, 210, replace=TRUE)
  
  
  ebu_forecast_1wk_nDA <- forecast_function(IC = rnorm(210,IC[,1], IC[,2]), ## sample IC
                                            FLARE = forecast_1wk,
                                            mu2 = parm_17$mu2,
                                            phi = parm_17$phi,
                                            omega = parm_17$omega,
                                            Q = parm_17$sd.pro,
                                            ndays = as.numeric(ndays))
  
  ebu_forecast_2wk_nDA <- forecast_function(IC = ebu_forecast_1wk_nDA,  ## sample IC
                                            FLARE = forecast_2wk,
                                            mu2 = parms$`pars[1]`,
                                            phi = parms$`pars[2]`,
                                            omega = parms$`pars[3]`,
                                            Q = parms$sd.pro,
                                            ndays = as.numeric(ndays))
  
  ebu_forecast_nDA <- rbind(ebu_forecast_1wk_nDA, ebu_forecast_2wk_nDA)%>% unname(.)%>%
    cbind(time,.)
  
  observed <- cbind(observed_date, t(rnorm(210,IC[,1], IC[,2])))
  
  ebu_forecast_nDA <- bind_rows(observed, ebu_forecast_nDA)
  
  ebu_forecast_raw_nDA <- melt(ebu_forecast_nDA, id = "time")
  
  ebu_forecast_nDA_summarized <- ebu_forecast_raw_nDA %>%
    group_by(time) %>%
    summarise(mean = mean(value),
              max = max(value),
              min = min(value),
              upper_95 = quantile(value, 0.95, na.rm = T),
              lower_95 = quantile(value, 0.05, na.rm = T),
              var = var(value))%>%
    mutate(forecast_date = start_forecast)
  saveRDS(ebu_forecast_nDA_summarized, paste0("./forecast_output/ebu_forecast_nDA_",dates[s],".rds"))
  
  
  
  jags.data.null = list(x_init = full_ebullition_model_alltrap_jags$ebu_rate[1],
                        Y = full_ebullition_model_alltrap_jags$ebu_rate, 
                        tau.obs = 1/((full_ebullition_model_alltrap_jags$ebu_rate_se)) ^ 2,
                        N = nrow(full_ebullition_model_alltrap_jags),
                        ndays = full_ebullition_model_alltrap_jags$ndays,
                        proc_prec = 1 / log(1 + (var(null_out_parms$sd.pro)/(mean(null_out_parms$sd.pro)^2))),
                        proc_mean = log((mean(null_out_parms$sd.pro)^2)/sqrt(mean(null_out_parms$sd.pro)^2 + var(null_out_parms$sd.pro))))
  
  nchain = 3
  chain_seeds <- c(200,800,1400)
  init <- list()
  for(i in 1:nchain){
    init[[i]] <- list(sd.pro = runif(1, 0.01,0.2),
                      .RNG.name = "base::Wichmann-Hill",
                      .RNG.seed = chain_seeds[i])
  }
  
  j.model   <- jags.model(file = null,
                          data = jags.data.null,
                          inits = init,
                          n.chains = 3)
  
  eval_ebu_null  <- coda.samples(model = j.model,
                                 variable.names = c("sd.pro"),
                                 n.iter = 200000, n.burnin = 20000, thin = 200)
  
  # Extract the parameter esimates from the ebullition model run for dates[s] and traps[h]
  ########################################################################################
  eval_ebu_null_parms <- eval_ebu_null %>%
    spread_draws(sd.pro) %>%
    rename(ensemble = .iteration) %>%
    mutate(forecast_date = start_forecast)%>%
    ungroup()%>%
    select(forecast_date, sd.pro)
  saveRDS(eval_ebu_null_parms, paste0("./forecast_output/null_ebullition_parameters_",dates[s],".rds"))
  
  time <- c(dates[s+1], dates[s+2])
  
  
  forecast_function <- function(IC, Q, ndays){
    est <- IC + rnorm(210, 0, sd = Q*ndays)
    return(est)
  }
  
  
  parms <- sample_n(eval_ebu_null_parms, 210, replace=TRUE)
  
  IC <- full_ebullition_model_alltrap_jags %>% filter(time == dates[s]) %>%
    select(ebu_rate, ebu_rate_se)
  
  
  observed_date <- full_ebullition_model_alltrap_jags %>% filter(time == dates[s]) %>%
    select(time)
  
  # initial condition uncertainty
  null_ebu_forecast_1wk <- forecast_function(IC = rnorm(210,IC[,1], IC[,2]), ## sample IC
                                             Q = parms$sd.pro,
                                             ndays = as.numeric(ndays))
  
  null_ebu_forecast_2wk <- forecast_function(IC = ebu_forecast_1wk,  ## sample IC
                                             Q = parms$sd.pro,
                                             ndays = as.numeric(ndays))
  
  null_ebu_forecast <- rbind(null_ebu_forecast_1wk, null_ebu_forecast_2wk)%>% unname(.)%>%
    cbind(time,.)
  
  
  observed <- cbind(observed_date, t(rnorm(210,IC[,1], IC[,2])))
  
  null_ebu_forecast <- bind_rows(observed, null_ebu_forecast)
  
  null_ebu_forecast_raw <- melt(null_ebu_forecast, id = "time")
  
  final_null_forecast_sumarized <- null_ebu_forecast_raw %>%
    group_by(time) %>%
    summarise(mean = mean(value),
              max = max(value),
              min = min(value),
              upper_95 = quantile(value, 0.95, na.rm = T),
              lower_95 = quantile(value, 0.05, na.rm = T),
              var = var(value))%>%
    mutate(forecast_date = start_forecast)
  
  saveRDS(final_null_forecast_sumarized, paste0("./forecast_output/null_ebu_forecast_",dates[s],".rds"))
  
}

loop <- proc.time()-t1
print(loop[3]/60)