# CH4cast
# Ryan McClure
# TRAINING MODELS WITH 2017 DATA

# 2017 MODEL TRAINING ----

#* TEMPERATURE SCALING MODEL ----
model.lm17 = ("lm17.txt")
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
                 file = model.lm17)

#* SS EBULLITION PREDICTION MODEL ----

model.ar17 = ("ar17_model.txt")
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
  }", file = model.ar17)

#* PERSISTENCE NULL MODEL ----
RandomWalk <- nimbleCode({
  
  #### Priors
  x[1] ~ dnorm(x_ic, sd = sd_ic)
  sd_add ~ dunif(0, 100)
  
  #### Process Model
  for(t in 2:n){
    pred[t] <- x[t-1]
    x[t] ~ dnorm(pred[t], sd = sd_add)
  }
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dnorm(x[t], sd = sd_obs)
  }
  
})

# RandomWalk = "
# model{
#   # Priors
#   x[1] ~ dnorm(x_ic,tau_init)
#   sd.pro ~ dunif(0, 1000)
#   tau.pro <-  pow(sd.pro, -2)
#   
#   #Informative priors on initial conditions based on first observation
#   predY[1] <- X[1]
#   Y[1] ~ dnorm(X[1], tau.obs[1])
#   
#   # Process Model
#   for(t in 2:n){
#     x[t] ~ dnorm(x[t-1],tau.pro)
#     Y[t] ~ dnorm(x[t],tau.obs[t])
#   }
# }"


# MCMC ANALYSES FOR 2017 ----

#* FILTER 2017 AND IMPUTE COVARIATE ----

full_ebullition_model_alltrap_jags <- full_ebullition_model_alltrap %>% 
  filter(time <= as.Date("2017-10-30"))%>%
  arrange(time)

for (i in colnames(full_ebullition_model_alltrap_jags[,c(2:10,12)])) {
  full_ebullition_model_alltrap_jags[,i] <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags[,i],option = "linear")   
}


constants <- list(n = length(full_ebullition_model_alltrap_jags$time),
                  x_ic = 0.0,
                  sd_ic = 0.1,
                  sd_obs = 1.0)

data <- list(y = full_ebullition_model_alltrap_jags$log_ebu_rate)

nchain = 3
inits <- list()
for(i in 1:nchain){
  y.samp = sample(full_ebullition_model_alltrap_jags$log_ebu_rate, length(full_ebullition_model_alltrap_jags$log_ebu_rate), replace = TRUE)
  inits[[i]] <- list(sd_add = sd(diff(y.samp)),
                     x = log(full_ebullition_model_alltrap_jags$log_ebu_rate))
}


nimble_out <- nimbleMCMC(code = RandomWalk,
                         data = data,
                         inits = inits,
                         constants = constants,
                         monitors = c("sd_add",
                                      "x",
                                      "y"),
                         niter = 11000,
                         nchains = 3,
                         samplesAsCodaMCMC = TRUE)







#* RUNJAGS FOR 2017 TEMP SCALE ----
jags.data.lm = list(Y= full_ebullition_model_alltrap_jags$hobo_temp,
                        N = nrow(full_ebullition_model_alltrap_jags),
                        C_mean = full_ebullition_model_alltrap_jags$water_temp_dam)

jags.params.lm.eval = c("sd.pro", "mu", "beta", "tau.pro")

j.lm.model   <- jags.model(file = model.lm17,
                           data = jags.data.lm,
                           n.chains = 3)

eval_temp  <- coda.samples(model = j.lm.model,
                      variable.names = jags.params.lm.eval,
                      n.iter = 10000, n.burnin = 1000)
plot(eval_temp)
print("TEMP SCALE MODEL DIAGNOSTICS")
print(gelman.diag(eval_temp))

#* RUNJAGS FOR 2017 EBULLITION MODEL ----

full_ebullition_model_alltrap_jags <- full_ebullition_model_alltrap_jags%>%
  select(time, log_ebu_rate, log_ebu_rate_sd, hobo_temp)%>%
  na.omit(.)

y_nogaps <- full_ebullition_model_alltrap_jags$log_ebu_rate[!is.na(full_ebullition_model_alltrap_jags$log_ebu_rate)]
full_ebullition_model_alltrap_jags$log_ebu_rate[is.nan(as.numeric(full_ebullition_model_alltrap_jags$log_ebu_rate))] <- 1.0

jags.data.ar = list(x_init = full_ebullition_model_alltrap_jags$log_ebu_rate[1],
                    Y = full_ebullition_model_alltrap_jags$log_ebu_rate, 
                    tau.obs = 1/((full_ebullition_model_alltrap_jags$log_ebu_rate_sd)/sqrt(4)) ^ 2,
                    N = nrow(full_ebullition_model_alltrap_jags), 
                    D = full_ebullition_model_alltrap_jags$hobo_temp)

nchain = 3
chain_seeds <- c(200,800,1400)
init <- list()
for(i in 1:nchain){
  init[[i]] <- list(sd.pro = 0.5,
                    mu2 = -5,
                    omega = 0.5, 
                    phi = 0.2,
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = chain_seeds[i])
}

j.model   <- jags.model(file = model.ar17,
                        data = jags.data.ar,
                        inits = init,
                        n.chains = 3)

eval_ebu  <- coda.samples(model = j.model,
                          variable.names = c("sd.pro", "mu2", "phi", "omega"),
                          n.iter = 30000, n.burnin = 10000)

plot(eval_ebu)
print("SS EBULLITION MODEL DIAGNOSTICS")
print(gelman.diag(eval_ebu))

#* RUNJAGS FOR 2017 PERSISTENCE NULL----
# Select site
site_data_var <- full_ebullition_model_alltrap_jags

#observed ebullition: Full time series with gaps
y_wgaps <- site_data_var$log_ebu_rate
sd_wgaps <- imputeTS::na_interpolation(site_data_var$log_ebu_rate_sd,option = "linear")
time <- c(site_data_var$time)
y_nogaps <- y_wgaps[!is.na(y_wgaps)]
y_wgaps_index <- 1:length(y_wgaps)
y_wgaps_index <- y_wgaps_index[!is.na(y_wgaps)]
init_x <- approx(x = time[!is.na(y_wgaps)], y = y_nogaps, xout = time, rule = 2)$y
data <- list(y = y_nogaps,
             y_wgaps_index = y_wgaps_index,
             nobs = length(y_wgaps_index),
             tau_obs = 1/(sd_wgaps ^ 2),
             n = length(y_wgaps),
             x_ic = 0.0)
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
j.model   <- jags.model (file = textConnection(RandomWalk),
                         data = data,
                         inits = init,
                         n.chains = 3)
jags.out   <- coda.samples(model = j.model,variable.names = c("tau_add","tau_init"), n.iter = 10000)
plot(jags.out)
print("PERSISTENCE NULL MODEL DIAGNOSTICS")
print(gelman.diag(jags.out))

# EXTRACT PARAMETER ESTIMATES ----

#* TEMPERATURE SCALING PARAMETERS ----
temp_out_parms <- eval_temp %>%
  spread_draws(sd.pro, mu, beta) %>%
  filter(.chain == 1) %>%
  rename(ensemble = .iteration) %>%
  ungroup()%>%
  select(sd.pro, mu, beta)%>%
  summarise(mean_sd.pro = mean(sd.pro),
            var_sd.pro = var(sd.pro),
            mean_mu = mean(mu),
            var.mu = var(mu),
            mean_beta = mean(beta),
            var_beta = var(beta))
print(temp_out_parms)

#* EBULLITION MODEL PARAMETERS ----
ebu_out_parms <- eval_ebu %>%
  spread_draws(sd.pro, mu2, phi, omega) %>%
  filter(.chain == 1) %>%
  rename(ensemble = .iteration) %>%
  ungroup()%>%
  select(sd.pro, mu2, phi, omega)%>%
  summarise(mean_sd.pro = mean(sd.pro),
            var_sd.pro =var(sd.pro),
            mean_mu2 = mean(mu2),
            var.mu2 = var(mu2),
            mean_omega = mean(omega),
            var_omega = var(omega),
            mean_phi = mean(phi),
            var_phi = var(phi))

print(ebu_out_parms)

#* PERSISTENCE NULL MODEL PARAMETERS ----
null_out_parms <- jags.out %>%
  spread_draws(tau_add, tau_init) %>%
  filter(.chain == 1) %>%
  rename(ensemble = .iteration) %>%
  summarise(mean_tau_add = mean(tau_add),
            var_tau_add = var(tau_add),
            mean_tau_init = mean(tau_init),
            var.tau_init = var(tau_init))

print(null_out_parms)