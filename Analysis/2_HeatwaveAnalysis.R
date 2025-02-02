rm(list = ls())
dev.off()

library(tidyverse)
library(heatwaveR)

setwd("D:/School/MichiganTech/Metabolism_Heatwave/Data")

wtemp_discharge_metab <- read_csv('wtemp_discharge_metab.csv')

# Run HW analysis ----

zz <- unique(wtemp_discharge_metab$site_no)

# Identify and summarize HW across all sites
for(i in 1:length(zz)){
  curDat = wtemp_discharge_metab[wtemp_discharge_metab$site_no == zz[i],]
  ts_Warm = ts2clm(curDat, x = date, y = Wtemp,
                   climatologyPeriod = c(min(curDat$date), max(curDat$date)))
  de_Warm = detect_event(ts_Warm, x = date, y = Wtemp )
  cat_Warm = category(de_Warm, y = Wtemp, S = FALSE)
  curEventsWarm = de_Warm$event
  curEventsWarm$Station = zz[i]
  curCatWarm = cat_Warm
  curCatWarm$Station = zz[i]
  if( i == 1){
    saveDatWarm = curEventsWarm
    saveCatWarm = curCatWarm
  } else{
    saveDatWarm = rbind(saveDatWarm, curEventsWarm)
    saveCatWarm = rbind(saveCatWarm, curCatWarm)
  }
}

NROW(saveDatWarm) # 2,129 events
round(mean(saveDatWarm$duration)) # 8 days
round(sd(saveDatWarm$duration)) # 5 days
max(saveDatWarm$duration) # 59 days
round(mean(saveDatWarm$intensity_max_relThresh),digits = 1) # 1.9 degrees C
round(max(saveDatWarm$intensity_max_relThresh),digits = 1) # 10.6 degrees C
round(mean(saveDatWarm$intensity_max),digits = 1) # 4.4 degrees C
round(max(saveDatWarm$intensity_max),digits = 1) # 13.5 degrees C

saveCatWarm %>%
  group_by(category) %>%
  summarise(Total_Events = n(), # Moderate = 1635 events, Strong = 460 events, Severe = 31 events, Extreme = 3 events
            Mean_Duration = round(mean(duration, na.rm = TRUE)), # Moderate = 8 days, Strong = 9 days, Severe = 12 days, Extreme = 10 days
            SD_Duration = round(sd(duration, na.rm = TRUE)),# Moderate = 5 days, Strong = 6 days, Severe = 7 days, Extreme = 5 days
            Max_Duration = max(duration, na.rm = TRUE))# Moderate = 53 days, Strong = 59 days, Severe = 37 days, Extreme = 15 days

# Extract water temperature and HW threshold during each HW event ----
for(i in 1:length(zz)){
  curDat = wtemp_discharge_metab[wtemp_discharge_metab$site_no == zz[i],]
  ts_Warm = ts2clm(curDat, x = date, y = Wtemp,
                   climatologyPeriod = c(min(curDat$date), max(curDat$date)))
  de_Warm = detect_event(ts_Warm, x = date, y = Wtemp)
  daily_intensity <- de_Warm$climatology
  heatwave_days <- daily_intensity[daily_intensity$event == TRUE, ]
  # daily_max_intensity <- heatwave_days[, c("t", "intensity")]
  curEventsWarm = heatwave_days
  curEventsWarm$Station = zz[i]
  if( i == 1){
    saveDatWarm2 = curEventsWarm
  } else{
    saveDatWarm2 = rbind(saveDatWarm2, curEventsWarm)
  }
}

saveDatWarm2 <- saveDatWarm2 %>%
  mutate(intensity_relThresh = Wtemp - thresh) %>%
  select(Station, date, Wtemp, seas, thresh, intensity_relThresh) %>%
  rename(site_no = Station)

saveDatWarm <- saveDatWarm %>%
  rename(site_no = Station)

saveCatWarm <- saveCatWarm %>%
  rename(site_no = Station)

# Data cleaning based on Appling et al. 2018
wtemp_discharge_metab <- wtemp_discharge_metab %>%
  mutate(GPP = ifelse(GPP < -0.5, NA, GPP),
         GPP = ifelse(GPP < 0, 0, GPP),
         ER = ifelse(ER > 0.5, NA, ER),
         ER = ifelse(ER > 0, 0, ER),
         NEP = GPP - abs(ER))

# Combine HW event metrics with HW categories
hw <- left_join(saveDatWarm,saveCatWarm, by = c('site_no','event_no')) %>%
  select(!c(event_name,peak_date,duration.y,index_start,index_peak,index_end)) %>%
  rename(duration = duration.x)

# Extract metabolism during heatwaves ----
subset_list <- list()

for (i in 1:nrow(hw)) {
  # Subset wtemp_discharge_metab based on the group and date range in the current row of hw
  subset_df <- wtemp_discharge_metab[wtemp_discharge_metab$site_no == hw$site_no[i] & wtemp_discharge_metab$date >= hw$date_start[i] & wtemp_discharge_metab$date <= hw$date_end[i], ]
  
  # Add the category column to the subset dataframe
  subset_df$category <- hw$category[i]
  subset_df$i_max <- hw$i_max[i]
  
  # Append the subset dataframe to the list
  subset_list[[i]] <- subset_df
}

# Combine all subsetted dataframes into a single dataframe
hw_metab <- do.call(rbind, subset_list)

hw_metab <- left_join(wtemp_discharge_metab,hw_metab, by = c('site_no','date')) %>%
  mutate(category = ifelse(is.na(category),'None',category),
         i_max = ifelse(is.na(i_max), 0, i_max)) %>%
  select(!c(long_name.y,lat.y,lon.y,Atemp.y,Wtemp.y,flow_cms.y,GPP.y,ER.y,NEP.y)) %>%
  rename(long_name = long_name.x,
         lat = lat.x,
         lon = lon.x,
         Atemp = Atemp.x,
         Wtemp = Wtemp.x,
         flow_cms = flow_cms.x,
         GPP = GPP.x,
         ER = ER.x,
         NEP = NEP.x) %>%
  mutate(category = case_match(category,
                               'None' ~ 'None',
                               'I Moderate' ~ 'Moderate',
                               'II Strong' ~ 'Strong',
                               'III Severe' ~ 'Severe',
                               'IV Extreme' ~ 'Extreme')) %>%
  as.data.frame()

hw_metab <- left_join(hw_metab, saveDatWarm2, by = c('site_no', 'date')) %>%
  select(!Wtemp.y) %>%
  rename(Wtemp = Wtemp.x) %>%
  mutate(intensity_relThresh = ifelse(is.na(intensity_relThresh), 0, intensity_relThresh),
         intensity_relThresh = ifelse(is.na(Wtemp), NA, intensity_relThresh),
         category = fct_relevel(category, 'None', 'Moderate', 'Strong', 'Severe', 'Extreme')) %>% 
  select(site_no, long_name, lat, lon, date, Wtemp, GPP, ER, NEP, category, i_max, seas, thresh, intensity_relThresh)

# What is the mean +/- sd & max HW intensity relative to the threshold for each severity classifications?
hw_metab %>%
  group_by(category) %>%
  summarise(Total_Days = n(),
            Mean_Intensity = round(mean(intensity_relThresh, na.rm = TRUE),2),
            SD_Intensity = round(sd(intensity_relThresh, na.rm = TRUE),2),
            Max_Intensity = round(max(intensity_relThresh, na.rm = TRUE),2))

# What is the mean +/- sd HW intensity relative to the threshold across all severity classifications?
hw_metab %>%
  filter(!category == 'None') %>%
  summarise(Mean_Intensity_All = round(mean(intensity_relThresh, na.rm = TRUE),2),
            SD_Intensity_All = round(sd(intensity_relThresh, na.rm = TRUE),2))

# What is the mean +/- sd & max HW intensity relative to the climatology for each severity classifications?
hw_metab %>%
  mutate(intensity_relSeas = Wtemp - seas,
         intensity_relSeas = ifelse(is.na(intensity_relSeas), 0, intensity_relSeas),
         intensity_relSeas = ifelse(is.na(Wtemp), NA, intensity_relSeas)) %>%
  group_by(category) %>%
  summarise(Total_Days = n(),
            Mean_Intensity_relSeas = round(mean(intensity_relSeas, na.rm = TRUE),2),
            SD_Intensity_relSeas = round(sd(intensity_relSeas, na.rm = TRUE),2),
            Max_Intensity_relSeas = round(max(intensity_relSeas, na.rm = TRUE),2))

# What is the mean +/- sd HW intensity relative to the climatology across all severity classifications?
hw_metab %>%
  mutate(intensity_relSeas = Wtemp - seas,
         intensity_relSeas = ifelse(is.na(intensity_relSeas), 0, intensity_relSeas),
         intensity_relSeas = ifelse(is.na(Wtemp), NA, intensity_relSeas)) %>%
  filter(!category == 'None') %>%
  summarise(Mean_Intensity_relSeas_All = round(mean(intensity_relSeas, na.rm = TRUE),2),
            SD_Intensity_relSeas_All = round(sd(intensity_relSeas, na.rm = TRUE),2))

write.csv(hw_metab, 'hw_metab.csv', row.names = FALSE)
