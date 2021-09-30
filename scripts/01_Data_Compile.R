#################################################################
# CH4castR                                                      #
# Ryan McClure                                                  #
# Data Wrangling Script                                         #
#################################################################


se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

####  FIRST TIME USERS ***NEED*** TO DOWNLOAD JAGS TO THEIR LOCAL COMPUTER #####

# Use command/control + click on the link to go directly to the download

# Jags download for Mac users can be found here: https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Mac%20OS%20X/
#* If you have a M1 arm 64 mac and are running the arm64 version of R 4.1.0. This will not work as JAGS is built on x86-64. 
#* For context, I learned this the hard way. 

# Jags download for PC users can be  found here: https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/

# Jags download for Linux users can be found here: https://packages.qa.debian.org/j/jags.html

# Download packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
pacman::p_load(tidyverse, MCMCvis, lubridate, tidybayes,
               ncdf4, reshape2, zoo, patchwork, hydroGOF, viridis,
               imputeTS, devtools, scales, forecast, coda, rjags, R2jags)

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
  mutate(ebu_rate = log(0.1 + Ebu_rate)) %>% 
  group_by(time) %>%
  summarize(ebu_rate_se = se(ebu_rate, na.rm = T),
            ebu_rate = mean(ebu_rate, na.rm = T)) %>%
  select(time, ebu_rate, ebu_rate_se)


ebu$ebu_rate[is.nan(as.numeric(ebu$ebu_rate))] <- NA

time <- as.data.frame(seq(from = as.Date("2017-05-07"), to = as.Date("2019-11-07"), by = "day"))%>%
  rename(time = `seq(from = as.Date(\"2017-05-07\"), to = as.Date(\"2019-11-07\"), by = \"day\")`)


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


