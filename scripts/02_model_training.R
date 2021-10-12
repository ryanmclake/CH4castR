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


#* RUNJAGS FOR 2017 TEMP SCALE ----


full_ebullition_model_alltrap_jags <- full_ebullition_model_alltrap%>%
   select(time, ebu_rate, ebu_rate_se, hobo_temp, water_temp_dam)%>%
   filter(time <= "2017-10-30")

for (i in colnames(full_ebullition_model_alltrap_jags[,c(3:5)])) {
   full_ebullition_model_alltrap_jags[,i] <- imputeTS::na_interpolation(full_ebullition_model_alltrap_jags[,i],option = "linear")   
}


full_ebullition_model_alltrap_jags <- na.omit(full_ebullition_model_alltrap_jags)%>%
   mutate(ndays = difftime(time,lag(time)))


nchain = 3
chain_seeds <- c(200,800,1400)
init <- list()
for(i in 1:nchain){
   init[[i]] <- list(sd.pro = runif(1, 0.5,1),
                     mu = runif(1,1,4),
                     beta = runif(1, 0.6, 0.8), 
                     .RNG.name = "base::Wichmann-Hill",
                     .RNG.seed = chain_seeds[i])
}

jags.data.lm = list(Y= full_ebullition_model_alltrap_jags$hobo_temp,
                    N = nrow(full_ebullition_model_alltrap_jags),
                    C_mean = full_ebullition_model_alltrap_jags$water_temp_dam)

jags.params.lm.eval = c("sd.pro", "mu", "beta")

j.lm.model   <- jags.model(file = model.lm17,
                           data = jags.data.lm,
                           n.chains = 3)

eval_temp  <- coda.samples(model = j.lm.model,
                           variable.names = jags.params.lm.eval,
                           n.iter = 200000, n.burnin = 20000, thin = 200)
#plot(eval_temp)
print("TEMP SCALE MODEL DIAGNOSTICS")
print(gelman.diag(eval_temp))


#* TEMPERATURE SCALING PARAMETERS ----

temp_out_parms <- eval_temp %>%
   spread_draws(sd.pro, mu, beta) %>%
   filter(.chain == 1)

print(temp_out_parms)


#* AR EBULLITION PREDICTION MODEL ----

model.ar17 = ("ar17_model.txt")
jagsscript = cat("
model {  
   
   #priors===================================================
   
   mu2 ~ dnorm(0,1e-6)
   sd.pro ~ dunif(0, 1000)
   phi ~ dnorm(0,1e-6)
   omega ~ dnorm(0,1e-6)
   
   #Informative priors on initial conditions based on first observation
   X[1] ~ dnorm(x_init, tau.obs[1])
   
   #end priors===============================================
   
   for(i in 2:N) {
      
      #process model=============================================
      
      tau.pro[i] <- 1/((sd.pro*ndays[i])*(sd.pro*ndays[i]))
      predX[i] <- mu2 + phi*X[i-1] + omega*D[i]
      X[i] ~ dnorm(predX[i],tau.pro[i])
      
      
      #end of process model======================================
      
      #data model================================================
      
      Y[i] ~ dnorm(X[i], tau.obs[i]) # Observation variation
            #end of data model=========================================
   }
  }", file = model.ar17)



#* RUNJAGS FOR 2017 EBULLITION MODEL ----


jags.data.ar = list(x_init = full_ebullition_model_alltrap_jags$ebu_rate[1],
                    Y = full_ebullition_model_alltrap_jags$ebu_rate, 
                    tau.obs = 1/((full_ebullition_model_alltrap_jags$ebu_rate_se)) ^ 2,
                    N = nrow(full_ebullition_model_alltrap_jags), 
                    D = full_ebullition_model_alltrap_jags$hobo_temp,
                    ndays = full_ebullition_model_alltrap_jags$ndays)

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

j.model   <- jags.model(file = model.ar17,
                        data = jags.data.ar,
                        inits = init,
                        n.chains = 3)

eval_ebu  <- coda.samples(model = j.model,
                          variable.names = c("sd.pro", "mu2", "phi", "omega"),
                          n.iter = 200000, n.burnin = 20000, thin = 200)

#plot(eval_ebu)
print("SS EBULLITION MODEL DIAGNOSTICS")
print(gelman.diag(eval_ebu))
plot(eval_ebu)


#* EBULLITION MODEL PARAMETERS ----

ebu_out_parms <- eval_ebu %>%
   spread_draws(sd.pro, mu2, phi, omega) %>%
   filter(.chain == 1)

print(ebu_out_parms)


############### HAS NOT BEEN UPDATED ######################
null = ("null_model.txt")
jagsscript = cat("
model {  
   
   #priors===================================================
   
   sd.pro ~ dunif(0, 1000)
   
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



#* RUNJAGS FOR 2017 EBULLITION MODEL ----


jags.data.null = list(x_init = full_ebullition_model_alltrap_jags$ebu_rate[1],
                      Y = full_ebullition_model_alltrap_jags$ebu_rate, 
                      tau.obs = 1/((full_ebullition_model_alltrap_jags$ebu_rate_se)) ^ 2,
                      N = nrow(full_ebullition_model_alltrap_jags),
                      ndays = full_ebullition_model_alltrap_jags$ndays)

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


plot(eval_ebu_null)
print("SS EBULLITION MODEL DIAGNOSTICS")
print(gelman.diag(eval_ebu_null))

#* PERSISTENCE NULL MODEL ----

#* PERSISTENCE NULL MODEL PARAMETERS ----
null_out_parms <- eval_ebu_null %>%
   spread_draws(sd.pro) %>%
   filter(.chain == 1)

print(null_out_parms)

