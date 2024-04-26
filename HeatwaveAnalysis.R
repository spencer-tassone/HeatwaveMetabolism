rm(list = ls())
dev.off()

library(tidyverse)
library(readr)
library(ggpubr)
library(ggsignif)
library(dataRetrieval)
library(heatwaveR)
# library(geodata)
# library(ggspatial)
# library(sf)
# install.packages("maptools", repos = "https://packagemanager.posit.co/cran/2023-10-13")
# library(remotes)
# install_github("JVAdams/jvamisc")
library(jvamisc)
library(patchwork)

# Metabolism Data Source
# https://www.sciencebase.gov/catalog/item/59bff507e4b091459a5e0982

setwd("D:/School/MichiganTech/Metabolism_Heatwave/")

metab <- read_tsv("daily_predictions.tsv")
site <- read_tsv("site_data.tsv")

metab <- metab %>%
  mutate(GPP = round(GPP,2),
         ER = round(ER, 2),
         temp.water = round(temp.water,2),
         discharge = round(discharge,2)) %>%
  select(date, site_name, GPP, ER) %>%
  left_join(., site, by = 'site_name') %>%
  select(date, nwis_id, GPP, ER) %>%
  rename(site_no = nwis_id)

# Download water temperature time series ----
# This download is time consuming, first run was saved and can be uploaded below
# stat_Cd <- "00003" # statistics parameter code = mean
# pcode_wtemp = '00010' # water temperature
# start_Date <- as.Date("1987-01-01")
# end_Date <- as.Date("2017-01-01")
# site_No <- unique(metab$site_no)
# 
# # Download daily mean water temperature time series ----
# wtemp <- readNWISdv(siteNumbers = site_No,
#                     parameterCd = pcode_wtemp,
#                     startDate = start_Date,
#                     endDate = end_Date,
#                     statCd = stat_Cd)%>%
#   renameNWISColumns()
# 
# wtemp <- wtemp %>%
#   mutate(nwis_id = paste('nwis_',site_no, sep = '')) %>%
#   select(c(nwis_id,site_no, Date, Wtemp, Wtemp_cd))

setwd("D:/School/MichiganTech/Metabolism_Heatwave")
# write.csv(wtemp,'wtemp.csv',row.names = F)
wtemp <- read.csv('wtemp.csv')
wtemp <- wtemp %>%
  select(!site_no) %>%
  separate(nwis_id,c('trash','site_no'),sep = '_') %>%
  select(!trash)

start_year <- wtemp %>%
  group_by(site_no) %>%
  summarise(start_year = min(year(Date)))

# Function to count datasets for each year
all_years <- data.frame(year = 1987:2017)

count_datasets <- function(year) {
  filtered_data <- start_year %>%
    filter(start_year <= year)  # Filter the dataset to include rows with start_date on or before the current year
  dataset_count <- nrow(filtered_data)  # Count the number of datasets
  return(dataset_count)
}

# Apply the count_datasets function to each year
result <- sapply(all_years$year, count_datasets)

# Create a data frame with the years and dataset counts
result_df <- data.frame(year = all_years$year, dataset_count = result)

result_df <- result_df %>%
  mutate(frac_avail = round(dataset_count/356,2),
         dataset_length = seq(30,0,-1),
         dataset_frac = round(dataset_length/30,2))

result_df %>%
  ggplot(aes(x = year, y = frac_avail)) +
  geom_point() +
  geom_point(aes(x = year, y = dataset_frac), col = 'red') +
  geom_vline(xintercept = 2002, linetype = 'longdash') +
  scale_x_continuous(breaks = seq(1987,2017,2)) +
  theme_bw()

wtemp_sites <- start_year %>%
  filter(start_year <= 1996)

Wtemp_daily <- wtemp %>%
  filter(site_no %in% wtemp_sites$site_no) %>%
  filter(!Date <= as.Date('1996-12-31')) # Daily mean water temperature for 20 years between 1997-2017

# Remove data that is not Approved (A), Approved Revised (A R), Approved Edited (A e) or Provisional (P) ----
Wtemp_daily$Wtemp[Wtemp_daily$Wtemp < 0] <- 0
table(Wtemp_daily$Wtemp_cd)
Wtemp_daily$Wtemp[Wtemp_daily$Wtemp_cd == "A [4]"] <- NA
Wtemp_daily$Wtemp[Wtemp_daily$Wtemp_cd == "A <"] <- NA

# Fill out time series for each site
startDate <- as.Date("1997-01-01")
endDate <- as.Date("2017-01-01")
x <- length(seq(from = startDate, to = endDate, by = 'day'))
y <- NROW(unique(Wtemp_daily$site_no)) # 57 sites
full_ts <- as.data.frame(rep(seq(from = startDate, to = endDate, by = "day"), times = y))
colnames(full_ts)[1] <- "Date"
full_site <- as.data.frame(rep(unique(Wtemp_daily$site_no),each = x))
colnames(full_site)[1] <- "site_no"
site_ts <- cbind(full_ts, full_site)
Wtemp_daily <- merge(site_ts, Wtemp_daily, by = c("site_no","Date"), all = TRUE)

# Linear interpolate for data gaps less than or equal to 2 days ----
sum(is.na(Wtemp_daily$Wtemp)) # 98,828

library(zoo)

Wtemp_daily <- Wtemp_daily %>%
  group_by(site_no) %>%
  mutate(Wtemp_int = na.approx(Wtemp, maxgap = 2, na.rm=F),
         Wtemp = ifelse(is.na(Wtemp), Wtemp_int, Wtemp),) %>%
  select(!Wtemp_int)

sum(is.na(Wtemp_daily$Wtemp)) # 97,188

# Model water temperature using meteorological data and multiple linear regression ---- 

site <- site %>%
  rename(site_no = nwis_id)

latlon <- left_join(Wtemp_daily,site, by = "site_no") %>%
  rename(site = site_no) %>%
  group_by(site) %>%
  summarise(lat = mean(lat),
            lon = mean(lon)) 

#* Grab meteorological data from Daymet web services to build water temperature multiple linear regressions ----
# https://daac.ornl.gov/
# https://www.nature.com/articles/s41597-021-00973-0#code-availability
library(daymetr)

write.table(latlon, paste0(tempdir(),"/latlon.csv"),
            sep = ",",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

met_dat <- download_daymet_batch(file_location = paste0(tempdir(),
                                                        "/latlon.csv"),
                                 start = 1997,
                                 end = 2017,
                                 internal = TRUE,
                                 simplify = TRUE)

met_dat_wide <- spread(met_dat, measurement, value)
z <- tail(which(met_dat_wide$site == '8181800'), 1)
met_dat_wide <- met_dat_wide %>%
  mutate(tmean = (tmax..deg.c. + tmin..deg.c.)/2,
         totalRadiation = (met_dat_wide$srad..W.m.2.*met_dat_wide$dayl..s.)/1000000, # calculation based on daymetr website https://daymet.ornl.gov/overview
         Date = as.Date(paste(year, yday, sep = "-"), "%Y-%j"),
         site_no = ifelse(row_number()<=z, paste0("0", site), site)) %>%
  filter(!Date >= as.Date('2017-01-02'))
met_dat_wide <- met_dat_wide[,c(18,17,3:4,8:16)]
# met_dat_wide$totalRadiation <- (met_dat_wide$srad..W.m.2.*met_dat_wide$dayl..s.)/1000000 
wmet <- merge(met_dat_wide, Wtemp_daily, by = c("site_no","Date"), all = T)

# remove leap days
remove_leap <- as.Date(c("1996-02-29","2000-02-29","2004-02-29",
                         "2008-02-29","2012-02-29","2016-02-29","2020-02-29"))
wmet <- wmet %>%
  filter(!Date %in% remove_leap)

# day of year that does not recognize leap day
wmet <- wmet %>% 
  mutate(DoY = day(Date),
         Month = month(Date),
         Year = year(Date)) %>% 
  group_by(Year, Month, site_no) %>%
  mutate(DoY = DoY - lag(DoY, default = 0)) %>%
  group_by(Year,site_no) %>%
  mutate(DoY = cumsum(DoY)) %>%
  select(-Month)

# Remove sites that have zero water temperature data
test <- wmet %>%
  group_by(site_no) %>%
  summarise(avail = all(is.na(Wtemp))) %>%
  filter(avail == TRUE) 

wmet <- wmet %>%
  filter(!site_no %in% test$site_no)

#* Fit models and remove sites with low predictability (R2 < 0.8) ----
library(broom)

models_fit <- wmet %>%
  group_by(site_no) %>%
  do(model = glance(lm(Wtemp~tmax..deg.c.+tmin..deg.c.+prcp..mm.day.+vp..Pa.+totalRadiation+DoY, data = .))) %>%
  unnest(model)

models_fit$r.squared <- round(models_fit$r.squared, digits = 2)
good_model_fits <- models_fit[models_fit$r.squared >= 0.80,] # 51 sites with r-square >= 0.80
round(mean(good_model_fits$r.squared),digits = 2) # answer is 0.92
round(sd(good_model_fits$r.squared),digits = 2) # answer is 0.04
wmet <- wmet %>%
  filter(site_no %in% good_model_fits$site_no)
met_dat_wide <- met_dat_wide %>%
  filter(site_no %in% good_model_fits$site_no)

#* Use models to get water temperature estimates ----
f <- function (.fit, .new_data) {
  predict(.fit, newdata = .new_data)
}

set.seed(8992)
wmet <- wmet %>%
  nest(data = -site_no) %>% 
  mutate(
    fit  = map(data, ~ lm(Wtemp~tmax..deg.c.+tmin..deg.c.+prcp..mm.day.+vp..Pa.+totalRadiation+DoY, data = .x)),
    yhat = map2(.x = fit, .y = data, f)
  ) %>% 
  unnest(cols = c(data, yhat)) %>% 
  select(-fit)

wmet <- wmet %>%
  mutate(yhat = ifelse(yhat < 0,0,round(yhat,2)),
         corWtemp = ifelse(is.na(Wtemp), round(yhat,2), Wtemp))

# Download daily mean discharge (Q) time series ----
pcode_discharge = '00060' # discharge

Q_daily_dat <- readNWISdv(siteNumbers = unique(wmet$site_no),
                          parameterCd = pcode_discharge,
                          startDate = startDate,
                          endDate = endDate) %>%
  renameNWISColumns()

# Remove tidal sites
tidal_sites <- Q_daily_dat %>%
  group_by(site_no) %>%
  filter(Flow < 0) %>%
  select(site_no) %>%
  distinct()

Q_daily_dat <- Q_daily_dat %>%
  filter(!site_no %in% tidal_sites$site_no)

wmet <- wmet %>%
  filter(!site_no %in% tidal_sites$site_no)
NROW(unique(wmet$site_no)) # 48 sites

# table(Q_daily_dat$Flow_cd)
# Q_daily_dat$Flow[Q_daily_dat$Flow_cd == "A e"] <- NA
Q_daily_dat <- Q_daily_dat %>%
  mutate(flow_cms = round(Flow*0.0283168,2)) %>%
  select(c(site_no,Date,flow_cms,Flow_cd))

# Determine how much Q data is missing for each station
yy <- NROW(unique(Q_daily_dat$site_no)) # 47 out of 48 sites have concurrent daily discharge data available
full_ts <- as.data.frame(rep(seq(from = startDate, to = endDate, by = "day"),times = yy)) 
colnames(full_ts)[1] <- "Date"
full_site <- as.data.frame(rep(unique(Q_daily_dat$site_no),each = x))
colnames(full_site)[1] <- "site_no"
site_ts <- cbind(full_ts, full_site)
Q_daily_dat <- merge(site_ts, Q_daily_dat, by = c("site_no","Date"), all = TRUE)

missing_data <- Q_daily_dat %>%
  group_by(site_no) %>%
  summarise(Total_Q_DataAvail = sum(!is.na(flow_cms)))
missing_data$Frac_Wtemp_Avail <- round(missing_data$Total_Q_DataAvail/x,2) # There are 9,497 days between 12/31/2021 - 1/1/1996

threshold <- 0.90
keep_sites_Q <- missing_data %>%
  filter(Frac_Wtemp_Avail >= threshold) # 40 out of 48 sites have enough discharge data

Q_daily_dat <- Q_daily_dat %>%
  filter(site_no %in% keep_sites_Q$site_no,
         !Date %in% remove_leap)

wtemp_discharge <- left_join(wmet,Q_daily_dat,by = c('site_no','Date'))

wtemp_discharge <- wtemp_discharge %>%
  mutate(nwis_id = site_no,
         date = Date) %>%
  select(c(nwis_id, date, tmean, corWtemp, flow_cms)) %>%
  rename(Atemp = tmean,
         Wtemp = corWtemp)
# wtemp_discharge <- wtemp_discharge[,c(1,20,2:3,11,18:19)]

wtemp_discharge_metab <- left_join(wtemp_discharge, metab, by = c('site_no', 'date'))

# Linear interpolate for metabolism & discharge data gaps less than or equal to 2 days ----
wtemp_discharge_metab <- wtemp_discharge_metab %>%
  group_by(site_no) %>%
  mutate(GPP_int = na.approx(GPP, maxgap = 2, na.rm = FALSE),
         ER_int = na.approx(ER, maxgap = 2, na.rm = FALSE),
         flow_cms_int = na.approx(flow_cms, maxgap = 2, na.rm = FALSE),
         GPP = ifelse(is.na(GPP), GPP_int, GPP),
         ER = ifelse(is.na(ER), ER_int, ER),
         flow_cms = ifelse(is.na(flow_cms), flow_cms_int, flow_cms)) %>%
  select(!c(GPP_int,ER_int,flow_cms_int))

sum(is.na(wtemp_discharge_metab$GPP)) # 234,224
sum(is.na(wtemp_discharge_metab$ER)) # 234,224
sum(is.na(wtemp_discharge_metab$flow_cms)) # 59,037

wtemp_discharge_metab <- left_join(wtemp_discharge_metab, site, by = "site_no") %>%
  select(site_no,long_name,lat,lon,date,Atemp,Wtemp,flow_cms,GPP,ER)

# lat_long <- wtemp_discharge_metab %>%
#   select(site_no, long_name, lat, lon) %>%
#   distinct()
# 
# write.csv(lat_long, '48USGS_site_locations.csv', row.names = FALSE)

# Run HW analysis ----

zz <- unique(wtemp_discharge_metab$site_no)
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
max(saveDatWarm$duration) # 59 days
round(mean(saveDatWarm$intensity_max_relThresh),digits = 1) # 1.9 degrees C
round(max(saveDatWarm$intensity_max_relThresh),digits = 1) # 10.6 degrees C
round(mean(saveDatWarm$intensity_max),digits = 1) # 4.4 degrees C
round(max(saveDatWarm$intensity_max),digits = 1) # 13.5 degrees C
table(saveCatWarm$category) # Moderate = 1635 events, Strong = 460 events, Severe = 31 events, Extreme = 3 events

# Mean GPP and ER during heatwaves

saveDatWarm <- saveDatWarm %>%
  rename(site_no = Station)

saveCatWarm <- saveCatWarm %>%
  rename(site_no = Station)

wtemp_discharge_metab$GPP[wtemp_discharge_metab$GPP < -0.5] <- NA
wtemp_discharge_metab$ER[wtemp_discharge_metab$ER > 0.5] <- NA
wtemp_discharge_metab$GPP[wtemp_discharge_metab$GPP < 0] <- 0
wtemp_discharge_metab$ER[wtemp_discharge_metab$ER > 0] <- 0

wtemp_discharge_metab <- wtemp_discharge_metab %>%
  mutate(NEM = GPP + ER)

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
  
  # Append the subset dataframe to the list
  subset_list[[i]] <- subset_df
}

# Combine all subsetted dataframes into a single dataframe
hw_metab <- do.call(rbind, subset_list)

hw_metab <- left_join(wtemp_discharge_metab,hw_metab, by = c('site_no','date')) %>%
  mutate(category = ifelse(is.na(category),'None',category)) %>%
  select(!c(long_name.y,lat.y,lon.y,Wtemp.y,flow_cms.y,GPP.y,ER.y,NEM.y)) %>%
  rename(long_name = long_name.x,
         lat = lat.x,
         lon = lon.x,
         Wtemp = Wtemp.x,
         flow_cms = flow_cms.x,
         GPP = GPP.x,
         ER = ER.x,
         NEM = NEM.x) %>%
  mutate(category = case_match(category,
                               'None' ~ 'None',
                               'I Moderate' ~ 'Moderate',
                               'II Strong' ~ 'Strong',
                               'III Severe' ~ 'Severe',
                               'IV Extreme' ~ 'Extreme')) %>%
  as.data.frame()

cols = c("None" = 'blue',
         "Moderate" = '#FFC866' ,
         "Strong" = '#FF6900' ,
         "Severe" = '#9E0000',
         "Extreme" = '#2D0000')

hw_metab$category <- factor(hw_metab$category,
                            levels = c('None',
                                       'Moderate',
                                       'Strong',
                                       'Severe',
                                       'Extreme'))

hw_metab %>%
  group_by(category) %>%
  summarise(TotalDays = n()) # None = 332,717 days, Moderate = 13,133 days, Strong = 4,209 days, Extreme = 31 days

# Figures (ECDF) ----
setwd("D:/School/MichiganTech/Metabolism_Heatwave")

gpp_plot <- ggplot(data = hw_metab, aes(x = GPP, color = category)) +
  stat_ecdf(geom = 'step', pad = F, linewidth = 1) +
  scale_color_manual(values = cols) +
  labs(color = 'Riverine\nHW Severity',
       x = expression(atop(Gross~Primary~Production,(g~O[2]~m^-2~d^-1))),
       y = 'Cumulative Distribution') +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  coord_cartesian(xlim=c(0, 15)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        legend.position = c(0.7,0.3),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

er_plot <- ggplot(data = hw_metab, aes(x = ER, color = category)) +
  stat_ecdf(geom = 'step', pad = F, linewidth = 1) +
  scale_color_manual(values = cols) +
  labs(x = expression(atop(Ecosystem~Respiration,(g~O[2]~m^-2~d^-1))),
       y = NULL) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  coord_cartesian(xlim=c(-20, 0)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        legend.position = 'none')

nem_plot <- ggplot(data = hw_metab, aes(x = NEM, color = category)) +
  stat_ecdf(geom = 'step', pad = F, linewidth = 1) +
  scale_color_manual(values = cols) +
  labs(x = expression(atop(Net~Ecosystem~Metabolism,(g~O[2]~m^-2~d^-1))),
       y = NULL) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  coord_cartesian(xlim=c(-16, 16)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        legend.position = 'none')

# width = 1500 height = 500
gpp_plot + er_plot + nem_plot + plot_layout(nrow = 1)

x <- hw_metab %>%
  filter(category == 'None') %>%
  summarise(Mean_GPP_None = round(mean(GPP, na.rm = T),2),
            Mean_ER_None = round(mean(ER, na.rm = T),2),
            Mean_NEM_None = round(mean(NEM, na.rm = T),2))

hw_metab %>%
  group_by(category) %>%
  summarise(nTotal = n(),
            Mean_GPP = mean(GPP, na.rm = T),
            SD_GPP = sd(GPP, na.rm = T),
            pct_diff_GPP = round(((Mean_GPP - x$Mean_GPP_None)/x$Mean_GPP_None)*100),
            Mean_ER = mean(ER, na.rm = T),
            SD_ER = sd(ER, na.rm = T),
            pct_diff_ER = round(((Mean_ER - x$Mean_ER_None)/x$Mean_ER_None)*100),
            Mean_NEM = mean(NEM, na.rm = T),
            SD_NEM = sd(NEM, na.rm = T),
            pct_diff_NEM = round(((Mean_NEM - x$Mean_NEM_None)/x$Mean_NEM_None)*100))

# HW Metabolism ANOVA ----
library(misty)

aov.b(GPP~category, data = hw_metab)
aov.b(ER~category, data = hw_metab)
aov.b(NEM~category, data = hw_metab)

# Figures (boxplot w/ sig) ----
gpp_sig <- hw_metab %>%
  ggplot(aes(x = forcats::fct_rev(category), y = GPP, color = category)) +
  geom_boxplot(outlier.shape = NA, aes(middle = mean(GPP))) +
  coord_flip(ylim=c(0, 15)) +
  scale_color_manual(values = cols) +
  annotate('segment', x = c(5,5,5,5), xend = c(4,3,2,1), y = c(11,12,13,14), yend = c(11,12,13,14), color = 'black') +
  annotate('text', x = c(4.5,4,3.5,3), y = c(11.5,12.5,13.5,14.5),
           label = c('< 0.001','NS','NS', 'NS'),
           color = 'black', angle = -90, size = 5) +
  labs(y = expression(atop(Gross~Primary~Production,(g~O[2]~m^-2~d^-1))),
       x = 'Heatwave Severity') +
  theme_bw() +
  theme(axis.title.x = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        legend.position = 'none',
        plot.margin = grid::unit(c(1,1,0,0), "mm"))

er_sig <- hw_metab %>%
  ggplot(aes(x = forcats::fct_rev(category), y = ER, color = category)) +
  geom_boxplot(outlier.shape = NA, aes(middle = mean(GPP))) +
  coord_flip(ylim=c(-20,6)) +
  scale_color_manual(values = cols) +
  annotate('segment', x = c(5,5,5,5), xend = c(4,3,2,1), y = c(0.5,2,3.5,5), yend = c(0.5,2,3.5,5), color = 'black') +
  annotate('text', x = c(4.5,4,3.5,3), y = c(1.25,2.75,4.25,5.75),
           label = c('< 0.001','< 0.001','NS', 'NS'),
           color = 'black', angle = -90, size = 5) +
  labs(y = expression(atop(Ecosystem~Respiration,(g~O[2]~m^-2~d^-1))),
       x = NULL) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 20, color = "black"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_blank(),
        legend.position = 'none',
        plot.margin = grid::unit(c(1,1,0,0), "mm"))

nem_sig <- hw_metab %>%
  ggplot(aes(x = forcats::fct_rev(category), y = NEM, color = category)) +
  geom_boxplot(outlier.shape = NA, aes(middle = mean(GPP))) +
  coord_flip(ylim=c(-16,15)) +
  scale_color_manual(values = cols) +
  annotate('segment', x = c(5,5,5,5), xend = c(4,3,2,1), y = c(8,9.75,11.5,13.25), yend = c(8,9.75,11.5,13.25), color = 'black') +
  annotate('text', x = c(4.5,4,3.5,3), y = c(9,10.75,12.5,14.25),
           label = c('NS','NS','< 0.001', '0.002'),
           color = 'black', angle = -90, size = 5) +
  labs(y = expression(atop(Net~Ecosystem~Metabolism,(g~O[2]~m^-2~d^-1))),
       x = NULL) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 20, color = "black"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_blank(),
        legend.position = 'none',
        plot.margin = grid::unit(c(1,1,0,0), "mm"))

# width = 1500 height = 500
gpp_sig + er_sig + nem_sig + plot_layout(nrow = 1)

# Q10 estimation ----
# width = 1000 height = 800
hw_metab %>%
  # filter(category == 'None') %>%
  mutate(Wtemp_round = round(Wtemp)) %>%
  group_by(Wtemp_round, category) %>%
  summarise(`Gross Primary Production` = median(GPP, na.rm = TRUE),
            `Ecosystem Respiration` = median(ER, na.rm = TRUE)) %>%
  pivot_longer(cols = c(`Gross Primary Production`, `Ecosystem Respiration`),
               names_to = 'MetabComponent',
               values_to = 'MedianMetabValue') %>%
  ggplot(aes(x = Wtemp_round, y = MedianMetabValue, color = category)) +
  geom_line(linewidth = 1) +
  labs(x = expression(Water~Temperature~(degree*C)),
       y = expression(Median~Metabolic~Rate~(g~O[2]~m^-2~d^-1)),
       color = 'Heatwave\nSeverity') +
  scale_color_manual(values = cols) +
  scale_y_continuous(breaks = seq(-10,15,5), limits = c(-12,18)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 14, colour = 'black'),
        legend.position = c(0.1,0.8),
        legend.text = element_text(size = 14, color = 'black'),
        legend.title = element_text(size = 14, color = 'black')) +
  facet_wrap(~MetabComponent)

q10_dat <- hw_metab %>%
  mutate(Wtemp_round = round(Wtemp)) %>%
  group_by(Wtemp_round, category) %>%
  summarise(MedianGPP = median(GPP, na.rm = TRUE),
            MedianER = median(ER, na.rm = TRUE)) %>%
  filter(Wtemp_round %in% c(10,30),
         category %in% c('None','Moderate','Strong'))

q10_dat

zz <- unique(q10_dat$category)
saveDat_Q10 <- NULL

for (i in zz) {
  curDat <- q10_dat[q10_dat$category == i,]
  
  MedianGPP10 <- curDat$MedianGPP[curDat$Wtemp_round == 10]
  MedianGPP30 <- curDat$MedianGPP[curDat$Wtemp_round == 30]
  MedianER10 <- curDat$MedianER[curDat$Wtemp_round == 10]
  MedianER30 <- curDat$MedianER[curDat$Wtemp_round == 30]
  Wtemp10 <- curDat$Wtemp_round[curDat$Wtemp_round == 10]
  Wtemp30 <- curDat$Wtemp_round[curDat$Wtemp_round == 30]
  
  q10_GPP <- round((MedianGPP30 / MedianGPP10) ^ (10 / (Wtemp30 - Wtemp10)),2)
  q10_ER <- round((MedianER30 / MedianER10) ^ (10 / (Wtemp30 - Wtemp10)),2)
  
  curInfo <- data.frame(category = i,
                        q10_GPP = q10_GPP,
                        q10_ER = q10_ER)
  
  if (is.null(saveDat)) {
    saveDat_Q10 <- curInfo
  } else {
    saveDat_Q10 <- rbind(saveDat_Q10, curInfo)
  }
}

saveDat_Q10
