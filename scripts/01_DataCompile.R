#################################################################
# CH4cast                                                       #
# Ryan McClure                                                  #
# Data Wrangling Script                                         #
#################################################################

if (!"pacman" %in% installed.packages()) install.packages("pacman")
pacman::p_load(tidyverse, rjags, runjags, MCMCvis, lubridate, tidybayes,
               R2jags, ncdf4, reshape2, zoo, patchwork, hydroGOF, viridis,
               imputeTS, devtools, lintr, scales)

### Pull together all of the observations ###

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Met data for the NOAA forecasts
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for more information on the meteorological data in EDI refer to this link:
# https://portal.edirepository.org/nis/mapbrowse?packageid=edi.271.5

### Specify the URL
inUrl1 <- "https://pasta.lternet.edu/package/data/eml/edi/389/5/3d1866fecfb8e17dc902c76436239431"

### Define as temporary file
infile1 <- tempfile()

### Download the file
try(download.file(inUrl1, infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1, infile1,method="auto")

### Read the downloaded file into R Environment
met <-read.csv(infile1,header=F,skip=1,sep=",",
               col.names=c("Site", 
                       "Reservoir",
                       "DateTime",
                       "Record",
                       "CR3000_Batt_V",
                       "CR3000Panel_temp_C",
                       "PAR_Average_umol_s_m2",
                       "PAR_Total_mmol_m2",     
                       "BP_Average_kPa",
                       "AirTemp_Average_C",
                       "RH_percent",
                       "Rain_Total_mm",
                       "WindSpeed_Average_m_s",
                       "WindDir_degrees",
                       "ShortwaveRadiationUp_Average_W_m2",
                       "ShortwaveRadiationDown_Average_W_m2",
                       "InfaredRadiationUp_Average_W_m2",
                       "InfaredRadiationDown_Average_W_m2",
                       "Albedo_Average_W_m2",
                       "Flag_PAR_Average_umol_s_m2",
                       "Note_PAR_Average_umol_s_m2",
                       "Flag_PAR_Total_mmol_m2",
                       "Note_PAR_Total_mmol_m2",
                       "Flag_BP_Average_kPa",
                       "Note_BP_Average_kPa",
                       "Flag_AirTemp_Average_C",
                       "Note_AirTemp_Average_C",
                       "Flag_RH_percent",
                       "Note_RH_percent",
                       "Flag_Rain_Total_mm",
                       "Note_Rain_Total_mm",
                       "Flag_WindSpeed_Average_m_s",
                       "Note_WindSpeed_Average_m_s",
                       "Flag_WindDir_degrees",
                       "Note_WindDir_degrees",
                       "Flag_ShortwaveRadiationUp_Average_W_m2",
                       "Note_ShortwaveRadiationUp_Average_W_m2",
                       "Flag_ShortwaveRadiationDown_Average_W_m2",
                       "Note_ShortwaveRadiationDown_Average_W_m2",
                       "Flag_InfaredRadiationUp_Average_W_m2",
                       "Note_InfaredRadiationUp_Average_W_m2",
                       "Flag_InfaredRadiationDown_Average_W_m2",
                       "Note_InfaredRadiationDown_Average_W_m2",
                       "Flag_Albedo_Average_W_m2",
                       "Note_Albedo_Average_W_m2"), 
               check.names=TRUE)

### Break the EDI link
unlink(infile1)

### Aggregate the met data to daily
met_sum <- met %>%
  select(DateTime, BP_Average_kPa, AirTemp_Average_C, WindSpeed_Average_m_s) %>%
  rename(time = DateTime) %>%
  mutate(time = as_date(time)) %>%
  group_by(time) %>%
  summarise(daily_bp = mean(BP_Average_kPa),
            daily_bp_sd = sd(BP_Average_kPa),
            daily_temp = mean(AirTemp_Average_C),
            daily_temp_sd = sd(AirTemp_Average_C),
            daily_wind = mean(WindSpeed_Average_m_s),
            daily_wind_sd = sd(WindSpeed_Average_m_s)) %>%
  mutate(diff_bp = daily_bp-lag(daily_bp)) %>%
  filter(time >= "2017-05-01") %>%
  filter(time <="2020-11-07") %>%
  select(time, daily_temp, daily_temp_sd, daily_wind, daily_wind_sd, diff_bp)

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
            mean_sd = sd(value))

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
            mean_sd = sd(temp))%>%
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
  dplyr::summarize(temperature_sd = sd(Sed_temp, na.rm = T),
                   temperature = mean(Sed_temp, na.rm = T))%>%
  select(time, temperature, temperature_sd)%>%
  filter(time > "2017-05-06")


# Read in the observed ebullition throughout the entire ebullition sampling period --> 2017 to 2019
# importantly, we are specifically selecting traps JUST FROM Transect 1
# Refer to Mcclure et al., 2020 for more information in regard to the other Transects in the reservoir

ebu <- read_csv("./observed/EDI_DATA_EBU_DIFF_DEPTH_2015_2019.csv") %>%
  rename(time = DateTime) %>%
  filter(Transect == "T1")%>%
  select(time, Site, Ebu_rate)%>%
  mutate(log_ebu_rate = log(Ebu_rate),
         log_ebu_rate = ifelse(is.infinite(log_ebu_rate),NA,log_ebu_rate))%>%
  group_by(time) %>%
  summarize(log_ebu_rate_sd = sd(log_ebu_rate, na.rm = T),
            log_ebu_rate = mean(log_ebu_rate, na.rm = T))%>%
  select(time, log_ebu_rate, log_ebu_rate_sd)


time <- as.data.frame(seq(from = as.Date("2017-05-07"), to = as.Date("2019-11-07"), by = "day"))%>%
  rename(time = `seq(from = as.Date(\"2017-05-07\"), to = as.Date(\"2019-11-07\"), by = \"day\")`)


# Join with all previous data and make the object that will feed directly into the Jags model
full_ebullition_model <- left_join(time, hobo, by = "time") %>%
  left_join(., met_sum, by = "time")%>%
  left_join(., water_temp, by = "time")%>%
  left_join(., ebu, by = "time")%>%
  mutate(year = year(time))%>%
  filter(year != "2018")
  
full_ebullition_model_17 <- full_ebullition_model%>%
  filter(time >= "2017-05-07") %>%
  filter(time <= "2017-10-30") %>%
  rename(water_temp_dam = mean) %>%
  rename(water_temp_dam_sd = mean_sd)

full_ebullition_model_19 <- full_ebullition_model%>%
  filter(time >= "2019-05-27") %>%
  filter(time <= "2019-11-07") %>%
  rename(water_temp_dam = mean) %>%
  rename(water_temp_dam_sd = mean_sd)

full_ebullition_model_alltrap <- bind_rows(full_ebullition_model_17, full_ebullition_model_19)%>%
  rename(hobo_temp = temperature, hobo_temp_sd = temperature_sd)

full_ebullition_model_alltrap$water_temp_dam[is.nan(as.numeric(full_ebullition_model_alltrap$water_temp_dam))] <- NA
full_ebullition_model_alltrap$water_temp_dam_sd[is.nan(as.numeric(full_ebullition_model_alltrap$water_temp_dam_sd))] <- NA
full_ebullition_model_alltrap$log_ebu_rate[is.nan(as.numeric(full_ebullition_model_alltrap$log_ebu_rate))] <- NA
full_ebullition_model_alltrap$log_ebu_rate_sd[is.nan(as.numeric(full_ebullition_model_alltrap$log_ebu_rate_sd))] <- NA
full_ebullition_model_alltrap$hobo_temp[is.nan(as.numeric(full_ebullition_model_alltrap$hobo_temp))] <- NA
full_ebullition_model_alltrap$hobo_temp_sd[is.nan(as.numeric(full_ebullition_model_alltrap$hobo_temp_sd))] <- NA

# Perform a check standard to confirm the Upstream Temps closely match the Hobo Temps
a <- ggplot(full_ebullition_model_alltrap, aes(water_temp_dam, hobo_temp))+
  geom_point()+
  geom_smooth(method = "lm")

