# testing nimble run
library(nimble)

ebu_raw <- read_csv("./observed/EDI_DATA_EBU_DIFF_DEPTH_2015_2019.csv") %>%
  rename(time = DateTime) %>%
  filter(Transect == "T1")%>%
  select(time, Site, Ebu_rate)%>%
  mutate(log_ebu_rate = log(Ebu_rate),
         log_ebu_rate = ifelse(is.infinite(log_ebu_rate),NA,log_ebu_rate))%>%
  select(-Ebu_rate)

timet1 <- as.data.frame(seq(from = as.Date("2017-05-07"), to = as.Date("2019-11-07"), by = "day"))%>%
  rename(time = `seq(from = as.Date(\"2017-05-07\"), to = as.Date(\"2019-11-07\"), by = \"day\")`)%>%
  mutate(Site = "T1e1")

timet2 <- as.data.frame(seq(from = as.Date("2017-05-07"), to = as.Date("2019-11-07"), by = "day"))%>%
  rename(time = `seq(from = as.Date(\"2017-05-07\"), to = as.Date(\"2019-11-07\"), by = \"day\")`)%>%
  mutate(Site = "T1e2")

timet3 <- as.data.frame(seq(from = as.Date("2017-05-07"), to = as.Date("2019-11-07"), by = "day"))%>%
  rename(time = `seq(from = as.Date(\"2017-05-07\"), to = as.Date(\"2019-11-07\"), by = \"day\")`)%>%
  mutate(Site = "T1e3")

timet4 <- as.data.frame(seq(from = as.Date("2017-05-07"), to = as.Date("2019-11-07"), by = "day"))%>%
  rename(time = `seq(from = as.Date(\"2017-05-07\"), to = as.Date(\"2019-11-07\"), by = \"day\")`)%>%
  mutate(Site = "T1e4")

water_tempt1 <- water_temp %>% mutate(Site = "T1e1")
water_tempt2 <- water_temp %>% mutate(Site = "T1e2")
water_tempt3 <- water_temp %>% mutate(Site = "T1e3")
water_tempt4 <- water_temp %>% mutate(Site = "T1e4")

hobot1 <- hobo %>% mutate(Site = "T1e1")
hobot2 <- hobo %>% mutate(Site = "T1e2")
hobot3 <- hobo %>% mutate(Site = "T1e3")
hobot4 <- hobo %>% mutate(Site = "T1e4")

water_temp_raw <- bind_rows(water_tempt1, water_tempt2, water_tempt3, water_tempt4)
time_raw <- bind_rows(timet1, timet2, timet3, timet4)
hobo_raw <- bind_rows(hobot1, hobot2, hobot3, hobot4)

full_ebullition_raw <- left_join(time_raw, water_temp_raw, by = c("time","Site")) %>%
  left_join(., hobo_raw, by = c("time","Site"))%>%
  left_join(., ebu_raw, by = c("time","Site"))%>%
  rename(hobo_temp = temperature, hobo_temp_sd = temperature_sd)%>%
  rename(water_temp_dam = mean) %>%
  rename(water_temp_dam_sd = mean_sd)%>%
  group_by(Site)%>%
  mutate(water_temp_dam = imputeTS::na_interpolation(water_temp_dam,option = "linear"),
         water_temp_dam_sd = imputeTS::na_interpolation(water_temp_dam_sd,option = "linear"),
         hobo_temp = imputeTS::na_interpolation(hobo_temp,option = "linear"),
         hobo_temp_sd = imputeTS::na_interpolation(hobo_temp_sd,option = "linear"))%>%
  arrange(time)%>%
  group_by(time, Site)%>%
  filter(row_number()==1)%>%
  ungroup()
  
group_factor <- factor(full_ebullition_raw$Site)
group_factor

group_index <- as.numeric(group_factor)
group_index

#Write Bayesian Models

l_m_equation_hierch <- nimbleCode({
  
  #Priors
  p1_global ~ dunif(0, 1)
  sigma_p1 ~ dunif(0.00001, 100)
  mu ~ dnorm(0,1)
  beta ~ dnorm(0,1)
  sd.pro ~ dunif(0.00001, 10000)
  
  for(j in 1:n_groups){
    p1[j] ~ dnorm(p1_global, sd = sigma_p1)
  }
  for(i in 1:n) {
    predX[i] <- mu + C[i]*beta; 
    X[i] ~ dnorm(predX[i], sd = sd.pro); # Process variation
    Y[i] ~ dnorm(X[i],  sd = sd.obs[i]); # Observation variation
    C[i] ~ dnorm(C_mean[i], sd = sd.pre[i]); # Covariate variation
  }
})


m_m_equation_hierch <- nimbleCode({
  
  #Priors
  p1_global ~ dunif(0, 20)
  sigma_p1 ~ dunif(0.00001, 100)
  
  sd.pro ~ dunif(0.00001, 100)
  mu2 ~ dnorm(0,5)
  phi ~ dnorm(0,1)
  omega ~ dnorm(0,5)
  
  for(j in 1:n_groups){
    p1[j] ~ dnorm(p1_global, sd = sigma_p1)
  }
    for(i in 2:n){
    predX[i] <- mu2 + phi*X[i-1] + omega*D[i]        ## process model
    X[i] ~ dnorm(predX[i], sd = sd.pro[i])          ## data model
    D[i] ~ dnorm(D_mean[i], sd = sd.obs[i])               ## Covariate model
  }
})


  
  # Select site and forecast date (#hashed out lines are for testing purposes)
    full_ebullition_raw_nimble <- full_ebullition_raw %>% 
    group_by(Site)%>%
    #filter(time <= dates[s])%>%
    filter(time <= as.Date("2019-07-01"))%>%
    arrange(time)

  # read the .nc file[h,s] from FLARE runs that correspond directly to dates ebullition from the traps were measured
  files <- list.files("./forecast_driver/")                                          
  start_forecast <- max(full_ebullition_raw_nimble$time)         
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


  # full_ebullition_model_alltrap_jags <- bind_rows(full_ebullition_model_alltrap_jags, forecast_jags[-1,])
  # full_ebullition_model_alltrap_jags$hobo_temp_sd <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags$hobo_temp_sd,option = "linear")
  # full_ebullition_model_alltrap_jags$log_ebu_rate_sd <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags$log_ebu_rate_sd,option = "linear")
  # 

  data = list(X = full_ebullition_raw_nimble$hobo_temp, 
                      sd.obs = full_ebullition_raw_nimble$hobo_temp_sd,
                      C_mean = full_ebullition_raw_nimble$water_temp_dam,
                      sd.pre = full_ebullition_raw_nimble$water_temp_dam_sd)
  

#Constants
  constants <- list(n = length(full_ebullition_raw_nimble$time),
                     n_groups = length(unique(group_index)),
                     group = group_index)

  ##Initial conditions for each chain
  nchain = 3
  inits <- list()
  for(i in 1:nchain){
    inits[[i]] <- list(p1_global = rnorm(1,3,0.1), 
                       mu = rnorm(1,0,0.1),
                       beta = rnorm(1,1,0.1),
                       sd.pro = runif(1, 0.2, 2),
                       sigma_p1 = runif(1, 0.2, 2))
  }
  
  #Run MCMC
  nimble.out <- nimbleMCMC(code = l_m_equation_hierch,
                           data = data,
                           inits = inits,
                           constants = constants,
                           monitors = c("p1_global", "mu", "beta", "sd.pro","sigma_p1", "p1", "Y"),
                           niter = 20000,
                           nchains = 3,
                           samplesAsCodaMCMC = TRUE)
  burnin <- 1000                               
  combined_fit <- window(nimble.out, start=burnin)
  gelman.diag(combined_fit)  ## determine convergence

  combined_fit %>%
    tidybayes::spread_draws(Y[Site]) %>%
    ggplot(aes(x = Y, y = factor(Site))) +
    stat_halfeye()  
  
  
  