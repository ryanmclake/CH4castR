# CH4castR user guide

<a href="url"><img src = "images/NSF.png" align="top" height="200" width="200" ></a>
<a href="url"><img src = "images/FLARE.jpg" align="top" height="200" width="200" ></a>
<a href="url"><img src = "images/CH4cast.png" align="top" height="200" width="180" ></a>

-----


:busts_in_silhouette: Ryan McClure, Quinn Thomas, Mary Lofton, Whitney Woelmer, and Cayelan Carey

:busts_in_silhouette: Special thanks to: Ashley Mickens, Bethany Bookout, Alexandria Hounshell, Abby Lewis, Heather Wander, Rachel Corrigan, James Maze, Dexter Howard, Jacob Wynne, and Nick Hammond, Paul Hanson, Erin Hotchkiss, Madeline Schreiber, and Bobbie Niederlehner

Questions?  :email: ryan333@vt.edu and rqthomas@vt.edu

## RUN CH4castR
Follow the myBinder link:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ryanmclake/CH4castR.git/HEAD)

-----

## Motivation

Thank you for checking out CH<sub>4</sub>castR. Near-term, ecological forecasting with iterative model refitting and uncertainty partitioning has great promise for improving our understanding of ecological processes and the predictive skill of ecological models, but to date has been infrequently applied to predict biogeochemical fluxes. Bubble fluxes of methane (CH<sub>4</sub>) from aquatic sediments to the atmosphere (ebullition) dominate freshwater greenhouse gas emissions, but it remains unknown how best to make robust near-term CH<sub>4</sub> ebullition predictions using models. Near-term forecasting workflows have the potential to address several current challenges in predicting CH<sub>4</sub> ebullition rates, including: development of models that can be applied across time horizons and ecosystems, identification of the timescales for which predictions can provide useful information, and quantification of uncertainty in predictions.

We developed and tested a near-term, iterative forecasting workflow of CH<sub>4</sub> ebullition rates in a small eutrophic reservoir throughout one open-water period. The workflow included the repeated updating of a CH<sub>4</sub> ebullition forecast model over time with newly-collected data via iterative model refitting. We compared the CH<sub>4</sub> forecasts from our workflow to both alternative forecasts generated without iterative model refitting and a persistence null model. This code provided is our workflow to generate forecasts of CH<sub>4</sub> ebullition rates in a eutrophic reservoir in Virginia, USA. 


## Prerequisites

CH<sub>4</sub>castR has been tested across Windows and Mac OS. 
*NOTE* --> CH<sub>4</sub>castR does not currently work with Apple's new M1 arm64 processor. Jags has not yet been prepped for this new processor. However, workarounds to get it working can be found at [SourceForge](https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/07e08a3605/)

### Jags will also need to be downloaded in order to run the forecasts
Jags download for Mac users can be found here: https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Mac%20OS%20X/
If you have a M1 arm 64 mac and are running the arm64 version of R 4.1.0, this code will not work as JAGS is built on x86-64 (I learned this the hard way!). 
Jags download for PC users can be found here: https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/
Jags download for Linux users can be found here: https://packages.qa.debian.org/j/jags.html

## Cloning NEON-forecast-code onto your computer (5 steps)
1. Go to the [CH<sub>4</sub>castR](https://github.com/ryanmclake/CH4castR) repository and copy the repo URL. 
2. Open R
3. Start a new project: File > New Project
4. Select: Version Control > Git
5. Paste the repo's URL into "Repository URL:", keep the project directory name as the default, select "open in new session", and click <b>New Project</b>

### Navigate to the "scritps" folder and open the following scripts
1. <i>01_Data_Compile.R</i>
2. <i>02_model_training.R</i>
3. <i>03_Generate_Forecasts.R</i>
4. <i>04_Forecast_Evaluations.R</i>

## Run Data Compiling Script
1. You can either source the <i>01_Data_Compile.R</i> script or run through it incrementally. 

### This script is setting up a few more configurations to execute a forecast. This includes:
Getting default packages
``` r
# Download packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
pacman::p_load(tidyverse, MCMCvis, lubridate, tidybayes,
               ncdf4, reshape2, zoo, patchwork, hydroGOF, viridis,
               imputeTS, devtools, scales, forecast, coda, rjags, R2jags)
```

Pull all of the available data from [Environmental Data Initiative](https://environmentaldatainitiative.org)
``` r
### Pull together all of the observations ###

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Catwalk data from dam site (See Figure 1 in MS)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for more information on the catwalk data in EDI refer to this link:
# https://portal.edirepository.org/nis/mapbrowse?packageid=edi.271.5

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/271/5/c1b1f16b8e3edbbff15444824b65fe8f"
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")
catwalk <-read.csv(infile1,header=F, skip=1, sep=",",
                   col.names=c("Reservoir",
                               "Site",
                               "DateTime",
                               "ThermistorTemp_C_surface",
                               "ThermistorTemp_C_1",
                               "ThermistorTemp_C_2",
                               "ThermistorTemp_C_3",
                               "ThermistorTemp_C_4",
                               "ThermistorTemp_C_5",
                               "ThermistorTemp_C_6",
                               "ThermistorTemp_C_7",
                               "ThermistorTemp_C_8",
                               "ThermistorTemp_C_9",
                               "RDO_mgL_5",
                               "RDOsat_percent_5",
                               "RDO_mgL_5_adjusted",
                               "RDOsat_percent_5_adjusted",
                               "RDOTemp_C_5",
                               "RDO_mgL_9",
                               "RDOsat_percent_9",
                               "RDO_mgL_9_adjusted",
                               "RDOsat_percent_9_adjusted",
                               "RDOTemp_C_9",
                               "EXOTemp_C_1",
                               "EXOCond_uScm_1",
                               "EXOSpCond_uScm_1",
                               "EXOTDS_mgL_1",
                               "EXODOsat_percent_1",
                               "EXODO_mgL_1",
                               "EXOChla_RFU_1",
                               "EXOChla_ugL_1",
                               "EXOBGAPC_RFU_1",
                               "EXOBGAPC_ugL_1",
                               "EXOfDOM_RFU_1",
                               "EXOfDOM_QSU_1",
                               "EXO_pressure",
                               "EXO_depth",
                               "EXO_battery",
                               "EXO_cablepower",
                               "EXO_wiper",
                               "Lvl_psi_9",
                               "LvlTemp_C_9",
                               "RECORD",
                               "CR6_Batt_V",
                               "CR6Panel_Temp_C",
                               "Flag_All",
                               "Flag_DO_1",
                               "Flag_DO_5",
                               "Flag_DO_9",
                               "Flag_Chla",
                               "Flag_Phyco",
                               "Flag_TDS",
                               "Flag_fDOM",
                               "Flag_Temp_Surf",
                               "Flag_Temp_1",
                               "Flag_Temp_2",
                               "Flag_Temp_3",
                               "Flag_Temp_4",
                               "Flag_Temp_5",
                               "Flag_Temp_6",
                               "Flag_Temp_7",
                               "Flag_Temp_8",
                               "Flag_Temp_9"),
                   check.names=TRUE)
unlink(infile1)

### Aggregate to daily temps between 2 and 3 meters
cat_sum <- catwalk %>%
  select(DateTime, ThermistorTemp_C_2, ThermistorTemp_C_3) %>%
  rename(time = DateTime)%>%
  mutate(time = as_date(time))%>%
  melt(., id = 'time')%>%
  mutate(value = as.numeric(value))%>%
  group_by(time)%>%
  summarise(mean = mean(value),
            mean_se = se(value))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### CTD data from dam site (See Figure 1 in MS)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for more information on the CTD data in EDI refer to this link:
# https://portal.edirepository.org/nis/mapbrowse?packageid=edi.200.11


inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/10/2461524a7da8f1906bfc3806d594f94c" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")

ctd <-read.csv(infile1, header=F, skip=1, sep=",",
               col.names=c("Reservoir", 
                           "Site",   
                           "Date",  
                           "Depth_m",
                           "Temp_C",
                           "DO_mgL",
                           "DO_pSat",
                           "Cond_uScm",
                           "Spec_Cond_uScm",
                           "Chla_ugL",
                           "Turb_NTU",     
                           "pH",
                           "ORP_mV",
                           "PAR_umolm2s",
                           "Desc_rate",
                           "Flag_Temp",
                           "Flag_DO",
                           "Flag_Cond",
                           "Flag_SpecCond",
                           "Flag_Chla",
                           "Flag_Turb",
                           "Flag_pH",
                           "Flag_ORP",
                           "Flag_PAR",
                           "Flag_DescRate"),
               check.names=TRUE)

unlink(infile1)

ctd_sum <- ctd %>% select(Reservoir, Date, Site, Depth_m, Temp_C)%>%
  filter(Reservoir == "FCR" & Site == 50)%>%
  filter(Date < "2018-07-05" & Date >= "2017-05-07")%>%
  filter(Depth_m > 2 & Depth_m <= 3)%>%
  rename(time = Date)%>%
  mutate(time = as_date(time))%>%
  select(time, Depth_m, Temp_C)%>%
  mutate(temp = as.numeric(Temp_C))%>%
  group_by(time)%>%
  summarise(mean = mean(temp),
            mean_se = se(temp))%>%
  filter(time != "2017-05-27" &    ### Extract CTD casts from days that FCR with Mixed using an epilimnetic Mixer.
         time != "2017-07-07" &    ### Refer to Lofton et al., 2018 and Chen et al., 2018 for details on the mixer
         time != "2017-07-08" & 
         time != "2017-07-09" &
         time != "2017-07-10" &
         time != "2017-07-11" &
         time != "2017-07-12" &
         time != "2017-10-08")

# make an observed damn site water temp column.
water_temp <- rbind(ctd_sum, cat_sum)%>%
  filter(time < "2019-11-08")

# Pull in hobo data 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hobo <- readRDS("./observed/EDI_DATA_HOBO_TEMPS_2017_2019.rds")%>%
  dplyr::rename(time = DateTime)%>%
  dplyr::mutate(time = as_date(time))%>%
  dplyr::group_by(time) %>%
  dplyr::summarize(temperature_se = se(Sed_temp, na.rm = T),
                   temperature = mean(Sed_temp, na.rm = T))%>%
  select(time, temperature, temperature_se)%>%
  filter(time > "2017-05-06")


# Read in the observed ebullition throughout the entire ebullition sampling period --> 2017 to 2019
# importantly, we are specifically selecting traps JUST FROM Transect 1
# Refer to Mcclure et al., 2020 for more information in regard to the other Transects in the reservoir

ebu <- read_csv("./observed/EDI_DATA_EBU_DIFF_DEPTH_2015_2019.csv") %>%
  rename(time = DateTime) %>%
  filter(Transect == "T1")%>%
  select(time, Site, Ebu_rate)%>%
  mutate(ebu_rate = ifelse(is.infinite(Ebu_rate),NA,Ebu_rate))%>%
  group_by(time) %>%
  summarize(ebu_rate_se = se(ebu_rate, na.rm = T),
            ebu_rate = mean(ebu_rate, na.rm = T))%>%
  select(time, ebu_rate, ebu_rate_se)

ebu$ebu_rate[is.nan(as.numeric(ebu$ebu_rate))] <- NA

time <- as.data.frame(seq(from = as.Date("2017-05-07"), to = as.Date("2019-11-07"), by = "day"))%>%
  rename(time = `seq(from = as.Date(\"2017-05-07\"), to = as.Date(\"2019-11-07\"), by = \"day\")`)
```

Join all of the data products that will feed into the JAGS model-fitting
``` r
# Join with all previous data and make the object that will feed directly into the Jags model
full_ebullition_model <- left_join(time, hobo, by = "time") %>%
  left_join(., water_temp, by = "time")%>%
  left_join(., ebu, by = "time")%>%
  mutate(year = year(time))%>%
  filter(year != "2018")
  
full_ebullition_model_17 <- full_ebullition_model%>%
  filter(time >= "2017-05-07") %>%
  filter(time <= "2017-10-30") %>%
  rename(water_temp_dam = mean) %>%
  rename(water_temp_dam_sd = mean_se)

full_ebullition_model_18 <- full_ebullition_model%>%
  filter(time >= "2018-05-07") %>%
  filter(time <= "2018-10-29") %>%
  rename(water_temp_dam = mean) %>%
  rename(water_temp_dam_sd = mean_se)

full_ebullition_model_19 <- full_ebullition_model%>%
  filter(time >= "2019-05-27") %>%
  filter(time <= "2019-11-07") %>%
  rename(water_temp_dam = mean) %>%
  rename(water_temp_dam_sd = mean_se)

full_ebullition_model_alltrap <- bind_rows(full_ebullition_model_17, full_ebullition_model_18, full_ebullition_model_19)%>%
  rename(hobo_temp = temperature, hobo_temp_se = temperature_se)

full_ebullition_model_alltrap$water_temp_dam[is.nan(as.numeric(full_ebullition_model_alltrap$water_temp_dam))] <- NA
full_ebullition_model_alltrap$water_temp_dam_sd[is.nan(as.numeric(full_ebullition_model_alltrap$water_temp_dam_sd))] <- NA
full_ebullition_model_alltrap$ebu_rate[is.nan(as.numeric(full_ebullition_model_alltrap$ebu_rate))] <- NA
full_ebullition_model_alltrap$ebu_rate_se[is.nan(as.numeric(full_ebullition_model_alltrap$ebu_rate_se))] <- NA
full_ebullition_model_alltrap$hobo_temp[is.nan(as.numeric(full_ebullition_model_alltrap$hobo_temp))] <- NA
full_ebullition_model_alltrap$hobo_temp_se[is.nan(as.numeric(full_ebullition_model_alltrap$hobo_temp_se))] <- NA
```

2. When script has finished running click on the <i>01_model_training.R</i> script. 

## Run Model Training Script
1. You can either source the <i>02_model_training.R</i> script or run through it incrementally. 

### This script is setting up the model-fitting for the 2019 forecast period
``` r
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
                    tau.obs = 1/((full_ebullition_model_alltrap_jags$ebu_rate_se)/sqrt(4)) ^ 2,
                    N = nrow(full_ebullition_model_alltrap_jags), 
                    D = full_ebullition_model_alltrap_jags$hobo_temp,
                    ndays = full_ebullition_model_alltrap_jags$ndays)

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
                          n.iter = 70000, n.burnin = 10000)

#plot(eval_ebu)
print("SS EBULLITION MODEL DIAGNOSTICS")
print(gelman.diag(eval_ebu))


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
                    tau.obs = 1/((full_ebullition_model_alltrap_jags$ebu_rate_se)/sqrt(4)) ^ 2,
                    N = nrow(full_ebullition_model_alltrap_jags),
                    ndays = full_ebullition_model_alltrap_jags$ndays)

nchain = 3
chain_seeds <- c(200,800,1400)
init <- list()
for(i in 1:nchain){
  init[[i]] <- list(sd.pro = 0.5,
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = chain_seeds[i])
}

j.model   <- jags.model(file = null,
                        data = jags.data.null,
                        inits = init,
                        n.chains = 3)

eval_ebu_null  <- coda.samples(model = j.model,
                          variable.names = c("sd.pro"),
                          n.iter = 70000, n.burnin = 10000)


plot(eval_ebu_null)
print("SS EBULLITION MODEL DIAGNOSTICS")
print(gelman.diag(eval_ebu_null))

#* PERSISTENCE NULL MODEL ----

#* PERSISTENCE NULL MODEL PARAMETERS ----
null_out_parms <- eval_ebu_null %>%
  spread_draws(sd.pro) %>%
  filter(.chain == 1)

print(null_out_parms)
```
2. When script has finished running click on the <i>03_Generate_Forecasts.R</i> script. 

## Run Forecasting Script

1. If everything has been properly aligned, you should also be able to click source and the script will run through. 
2. NOTE --> You might get an error at the end when sourcing. This is OK. It is still developing all the forecasts needed for the output.
3. This script includes the bulk for CH<sub>4</sub>castR. Everything else so far has been just staging to prep for the foreacsts. The coding block is almost the same as the model training script above. 
4. The forecasts are all saved in the forecast_output folder as .rds files. 

## Run Forecast Evaluation Script

1. If everything has been properly aligned, you should also be able to click source and the script will run through. 
2. This will generate plots in the figures folder. This includes Figure 3, 4, 5, and 6 in the MS. Figures 1 and 2 have been manually added to the figures folder for transparency and access. 


# Thanks again and Happy Forecasting!
