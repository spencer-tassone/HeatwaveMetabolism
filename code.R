library(readr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(lubridate)
library(dataRetrieval)
library(heatwaveR)
library(geodata)
library(ggspatial)
library(sf)
# library(remotes)
# install_github("JVAdams/jvamisc")
library(jvamisc)

rm(list = ls())
dev.off()

# Metabolism Data Source
# https://www.sciencebase.gov/catalog/item/59bff507e4b091459a5e0982

metab <- read_tsv("D:/School/MichiganTech/Metabolism_Heatwave/daily_predictions.tsv")
site <- read_tsv("D:/School/MichiganTech/Metabolism_Heatwave/site_data.tsv")

metab <- metab %>%
  mutate(GPP = round(GPP,2),
         ER = round(ER, 2),
         temp.water = round(temp.water,2),
         discharge = round(discharge,2)) %>%
  select(date,site_name,GPP,ER) %>%
  left_join(.,site,by = 'site_name') %>%
  select(date,nwis_id,GPP,ER)
  
stat_Cd <- "00003" # statistics parameter code = mean
pcode_wtemp = '00010' # water temperature
pcode_discharge = '00060' # discharge
start_Date <- as.Date("1987-01-01")
end_Date <- as.Date("2017-01-01")
site_No <- unique(metab$nwis_id)

# Download daily mean water temperature time series ----
wtemp <- readNWISdv(siteNumbers = site_No,
                    parameterCd = pcode_wtemp,
                    startDate = start_Date,
                    endDate = end_Date,
                    statCd = stat_Cd)
wtemp <- renameNWISColumns(wtemp)
wtemp <- wtemp %>%
  select(c(site_no, Date, Wtemp, Wtemp_cd))

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
         dataset_length = seq(31,1,-1),
         dataset_frac = round(dataset_length/31,2))

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
  filter(!Date <= as.Date('1995-12-31'))

# Remove data that is not Approved (A), Approved Revised (A R), Approved Edited (A e) or Provisional (P) ----
Wtemp_daily$Wtemp[Wtemp_daily$Wtemp >= 50] <- NA
Wtemp_daily$Wtemp[Wtemp_daily$Wtemp < 0] <- 0
Wtemp_daily$Wtemp[Wtemp_daily$Wtemp_cd == "A [4]"] <- NA
Wtemp_daily$Wtemp[Wtemp_daily$Wtemp_cd == "A <"] <- NA
Wtemp_daily$Wtemp[Wtemp_daily$Wtemp_cd == "P ***"] <- NA
Wtemp_daily$Wtemp[Wtemp_daily$Wtemp_cd == "P [4]"] <- NA
Wtemp_daily$Wtemp[Wtemp_daily$Wtemp_cd == "P Dis"] <- NA
Wtemp_daily$Wtemp[Wtemp_daily$Wtemp_cd == "P Eqp"] <- NA
Wtemp_daily$Wtemp[Wtemp_daily$Wtemp_cd == "P Mnt"] <- NA

# Fill out time series for each site
startDate <- as.Date("1996-01-01")
endDate <- as.Date("2017-01-01")
NROW(unique(Wtemp_daily$site_no)) # 57
full_ts <- as.data.frame(rep(seq(from = startDate, to = endDate, by = "day"), times = 57))
colnames(full_ts)[1] <- "Date"
length(seq(from = startDate, to = endDate, by = 'day'))
full_site <- as.data.frame(rep(unique(Wtemp_daily$site_no),each = 7672))
colnames(full_site)[1] <- "site_no"
site_ts <- cbind(full_ts, full_site)
Wtemp_daily <- merge(site_ts, Wtemp_daily, by = c("site_no","Date"), all = TRUE)

# Linear interpolate for data gaps less than or equal to 2 days
sum(is.na(Wtemp_daily$Wtemp)) # 107,986

library(zoo)

Wtemp_daily <- Wtemp_daily %>%
  group_by(site_no) %>%
  mutate(Wtemp_int = na.approx(Wtemp, maxgap = 2, na.rm=F),
         Wtemp = ifelse(is.na(Wtemp), Wtemp_int, Wtemp),) %>%
  select(!Wtemp_int)

sum(is.na(Wtemp_daily$Wtemp)) # 106,278

detach("package:zoo", unload = TRUE)

# Determine how much Wtemp data is missing for each station
missing_data <- Wtemp_daily %>%
  group_by(site_no) %>%
  summarise(Total_Wtemp_DataAvail = sum(!is.na(Wtemp)))

missing_data$Frac_Wtemp_Avail <- round(missing_data$Total_Wtemp_DataAvail/7672,2) # There are 7,672 days between 2017-01-01 - 1996-01-01

threshold <- 0.85
keep_sites_temp <- missing_data %>%
  filter(Frac_Wtemp_Avail >= threshold) # 30 stations @ 85% data availability

Wtemp_daily <- Wtemp_daily %>%
  filter(site_no %in% keep_sites_temp$site_no)

# Download daily mean discharge (Q) time series ----
start_Date <- as.Date("1996-01-01")
end_Date <- as.Date("2017-01-01")

Q_daily_dat <- readNWISdv(siteNumbers = keep_sites_temp$site_no,
                          parameterCd = pcode_discharge,
                          startDate = startDate,
                          endDate = endDate)
Q_daily_dat <- renameNWISColumns(Q_daily_dat)

# Remove tidal sites
tidal_sites <- Q_daily_dat %>%
  group_by(site_no) %>%
  filter(Flow < 0) %>%
  select(site_no) %>%
  distinct()

Q_daily_dat <- Q_daily_dat %>%
  filter(!site_no %in% tidal_sites$site_no)

Wtemp_daily <- Wtemp_daily %>%
  filter(!site_no %in% tidal_sites$site_no)

# table(Q_daily_dat$Flow_cd)
# Q_daily_dat$Flow[Q_daily_dat$Flow_cd == "A e"] <- NA
Q_daily_dat <- Q_daily_dat %>%
  mutate(flow_cms = round(Flow*0.0283168,2)) %>%
  select(c(site_no,Date,flow_cms,Flow_cd))

# Determine how much Q data is missing for each station
startDate <- as.Date("1996-01-01")
endDate <- as.Date("2017-01-01")
NROW(unique(Q_daily_dat$site_no)) # 29 out of 29 sites have concurrent daily discharge data available
full_ts <- as.data.frame(rep(seq(from = startDate, to = endDate, by = "day"),times = 29)) 
colnames(full_ts)[1] <- "Date"
full_site <- as.data.frame(rep(unique(Q_daily_dat$site_no),each = 7672))
colnames(full_site)[1] <- "site_no"
site_ts <- cbind(full_ts, full_site)
Q_daily_dat <- merge(site_ts, Q_daily_dat, by = c("site_no","Date"), all = TRUE)

missing_data <- Q_daily_dat %>%
  group_by(site_no) %>%
  summarise(Total_Q_DataAvail = sum(!is.na(flow_cms)))
missing_data$Frac_Wtemp_Avail <- round(missing_data$Total_Q_DataAvail/7672,2) # There are 9,497 days between 12/31/2021 - 1/1/1996

threshold <- 0.85
keep_sites_Q <- missing_data %>%
  filter(Frac_Wtemp_Avail >= threshold) # 27 out of 29 sites have enough discharge data

Q_daily_dat <- Q_daily_dat %>%
  filter(site_no %in% keep_sites_Q$site_no)

wtemp_discharge <- left_join(Wtemp_daily,Q_daily_dat,by = c('site_no','Date'))
wtemp_discharge <- wtemp_discharge %>%
  mutate(nwis_id = site_no,
         date = Date) %>%
  select(!c(Flow_cd,Wtemp_cd,site_no,Date))
wtemp_discharge <- wtemp_discharge[,c(4:5,2:3)]

wtemp_discharge_metab <- left_join(wtemp_discharge, metab, by = c('date','nwis_id'))

# Linear interpolate for data gaps less than or equal to 2 days
sum(is.na(wtemp_discharge_metab$GPP)) # 142,097
sum(is.na(wtemp_discharge_metab$ER)) # 142,097
sum(is.na(wtemp_discharge_metab$flow_cms)) # 17,547

library(zoo)

wtemp_discharge_metab <- wtemp_discharge_metab %>%
  group_by(nwis_id) %>%
  mutate(GPP_int = na.approx(GPP, maxgap = 2, na.rm=F),
         ER_int = na.approx(ER, maxgap = 2, na.rm=F),
         flow_cms_int = na.approx(flow_cms, maxgap = 2, na.rm=F),
         GPP = ifelse(is.na(GPP), GPP_int, GPP),
         ER = ifelse(is.na(ER), ER_int, ER),
         flow_cms = ifelse(is.na(flow_cms), flow_cms_int, flow_cms)) %>%
  select(!c(GPP_int,ER_int,flow_cms_int))

sum(is.na(wtemp_discharge_metab$GPP)) # 139,041
sum(is.na(wtemp_discharge_metab$ER)) # 139,041
sum(is.na(wtemp_discharge_metab$flow_cms)) # 17,545

detach("package:zoo", unload = TRUE)

wtemp_discharge_metab <- left_join(wtemp_discharge_metab, site, by = "nwis_id") %>%
  select(nwis_id,long_name,lat,lon,date,Wtemp,flow_cms,GPP,ER)

# Map ----

lat_long <- wtemp_discharge_metab %>%
  rename(site_no = nwis_id) %>%
  group_by(site_no) %>%
  summarise(lat = mean(lat),
            lon = mean(lon)) %>%
  mutate(QyesNo = ifelse(site_no %in% keep_sites_Q$site_no, 'Yes', 'No'))

sf_all <- st_as_sf(lat_long, 
                   coords = c("lon", "lat"),
                   crs = 4269)

can <- gadm(country = 'CAN', level = 1, path=tempdir()) %>%
  st_as_sf()
usa <- gadm(country = 'USA', level = 1, path=tempdir()) %>%
  st_as_sf()
mex <- gadm(country = 'MEX', level = 1, path=tempdir()) %>%
  st_as_sf()

lakes <- rnaturalearth::ne_download(scale = 110, 
                                    type = 'lakes', 
                                    category = 'physical') %>% 
  st_as_sf(lakes110, crs = 4269) # crs = 4326 (WGS84), crs = 4269(NAD83)

cols <- c("Yes" = "dodgerblue", "No" = "firebrick2")

ggplot() +
  geom_sf(data = can,
          mapping = aes(geometry = geometry),
          color = "black",
          fill = 'white') +
  geom_sf(data = usa,
          mapping = aes(geometry = geometry),
          color = "black",
          fill = 'white') +
  geom_sf(data = mex,
          mapping = aes(geometry = geometry),
          color = "black",
          fill = 'white') +
  geom_sf(data = lakes,
          mapping = aes(geometry = geometry),
          color = "black",
          fill = "white")  +
  geom_sf(data = sf_all,
          aes(fill = QyesNo),
          color = "black",
          shape = 21,
          size = 3,
          alpha = 0.6) +
  scale_fill_manual(labels = c("No (n = 2)",
                               "Yes (n = 27)"),
                    values = cols) +
  labs(fill = "Discharge Available") +
  coord_sf(xlim = c(-125.5, -66), ylim = c(24.3, 49.5), expand = F) +
  annotation_scale(location = "br",
                   width_hint = 0.1, 
                   text_cex = 1) +
  annotation_north_arrow(location = "br", 
                         which_north = "grid", 
                         pad_x = unit(0.15, "in"),
                         pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        legend.position = c(0.9,0.3),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.title.align = 0.5,
        legend.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"))

lat_long <- as.data.frame(lat_long[,c(3,2,1,4)])
lat_long$state <- latlong2(lat_long, to = 'state')
lat_long$state <- stringr::str_to_title(lat_long$state)
lat_long$state_abb <- state.abb[match(lat_long$state, state.name)]
# View(lat_long)

# Run HW analysis ----

zz <- unique(wtemp_discharge_metab$nwis_id)
for(i in 1:length(zz)){
  curDat = wtemp_discharge_metab[wtemp_discharge_metab$nwis_id == zz[i],]
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

NROW(saveDatWarm) # 1,396 events
round(mean(saveDatWarm$duration)) # 9 days
max(saveDatWarm$duration) # 131 days
round(mean(saveDatWarm$intensity_max_relThresh),digits = 1) # 1.8 degrees C
round(max(saveDatWarm$intensity_max_relThresh),digits = 1) # 9.5 degrees C
round(mean(saveDatWarm$intensity_max),digits = 1) # 4.0 degrees C
round(max(saveDatWarm$intensity_max),digits = 1) # 12.6 degrees C

hw_site <- saveDatWarm %>%
  group_by(Station) %>%
  summarise(Total = max(event_no))

# Mean GPP and ER during heatwaves

saveDatWarm <- saveDatWarm %>%
  rename(site_no = Station)

wtemp_discharge_metab <- wtemp_discharge_metab %>%
  rename(site_no = nwis_id)

sum(wtemp_discharge_metab$GPP < 0 & wtemp_discharge_metab$GPP >= -0.5, na.rm = T) # 6,108
sum(wtemp_discharge_metab$ER > 0 & wtemp_discharge_metab$ER <= 0.5, na.rm = T) # 728

wtemp_discharge_metab$GPP[wtemp_discharge_metab$GPP < -0.5] <- NA
wtemp_discharge_metab$ER[wtemp_discharge_metab$ER > 0.5] <- NA
wtemp_discharge_metab$GPP[wtemp_discharge_metab$GPP < 0] <- 0
wtemp_discharge_metab$ER[wtemp_discharge_metab$ER > 0] <- 0

wtemp_discharge_metab <- wtemp_discharge_metab %>%
  mutate(NEM = GPP + ER)

ranges <- mapply(function(x, y, z)  seq.Date(y, z, 1), saveDatWarm$site_no,  saveDatWarm$date_start, saveDatWarm$date_end, USE.NAMES = TRUE)
saveDatWarm$Mean_GPP <- mapply(function(a, b)
  mean(wtemp_discharge_metab$GPP[wtemp_discharge_metab$site_no == b][match(a, wtemp_discharge_metab$date[wtemp_discharge_metab$site_no == b])], na.rm = T), ranges, names(ranges))
saveDatWarm$Mean_ER <- mapply(function(a, b)
  mean(wtemp_discharge_metab$ER[wtemp_discharge_metab$site_no == b][match(a, wtemp_discharge_metab$date[wtemp_discharge_metab$site_no == b])], na.rm = T), ranges, names(ranges))
saveDatWarm$Mean_NEM <- mapply(function(a, b)
  mean(wtemp_discharge_metab$NEM[wtemp_discharge_metab$site_no == b][match(a, wtemp_discharge_metab$date[wtemp_discharge_metab$site_no == b])], na.rm = T), ranges, names(ranges))

saveCatWarm <- saveCatWarm %>%
  rename(site_no = Station)

hw <- left_join(saveDatWarm,saveCatWarm, by = c('site_no','event_no')) %>%
  select(!c(event_name,peak_date,duration.y,index_start,index_peak,index_end)) %>%
  rename(duration = duration.x)

# Extract metabolism during heatwaves ----
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
  as.data.frame()

cols = c("None" = 'blue',
         "I Moderate" = '#FFC866' ,
         "II Strong" = '#FF6900' ,
         "III Severe" = '#9E0000',
         "IV Extreme" = '#2D0000')

hw_metab$category <- factor(hw_metab$category,
                            levels = c('None',
                                       'I Moderate',
                                       'II Strong',
                                       'III Severe',
                                       'IV Extreme'))

# Figures (ECDF) ----
setwd("D:/School/MichiganTech/Metabolism_Heatwave")

# width = 800 height = 600
nem_plot <- ggplot(data = hw_metab, aes(x = NEM, color = category)) +
  stat_ecdf(geom = 'step', pad = F, size = 1) +
  scale_color_manual(values = cols) +
  labs(color = 'Riverine\nHW Severity') +
  xlab(bquote('Net Ecosystem Metabolism (g'~O[2]~ m^-2~d^-1*')')) +
  ylab('ECDF') +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        legend.position = 'none')
        # legend.position = c(0.15,0.8),
        # legend.title = element_text(size = 16),
        # legend.text = element_text(size = 16))
nem_plot

# width = 800 height = 600
er_plot <- ggplot(data = hw_metab, aes(x = ER, color = category)) +
  stat_ecdf(geom = 'step', pad = F, size = 1) +
  scale_color_manual(values = cols) +
  labs(color = 'Riverine\nHW Severity') +
  xlab(bquote('Ecosystem Respiration (g'~O[2]~ m^-2~d^-1*')')) +
  ylab('ECDF') +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        legend.position = 'none')
        # legend.position = c(0.15,0.8),
        # legend.title = element_text(size = 16),
        # legend.text = element_text(size = 16))
er_plot

# width = 800 height = 600
gpp_plot <- ggplot(data = hw_metab, aes(x = GPP, color = category)) +
  stat_ecdf(geom = 'step', pad = F, size = 1) +
  scale_color_manual(values = cols) +
  labs(color = 'Riverine\nHW Severity') +
  xlab(bquote('Gross Primary Production (g'~O[2]~ m^-2~d^-1*')')) +
  ylab('ECDF') +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        legend.position = c(0.85,0.35),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))
gpp_plot

# width = 600 height = 1200
ggarrange(gpp_plot,er_plot,nem_plot, ncol = 1)

gpp_plot2 <- ggplot(data = hw_metab, aes(x = GPP, color = category)) +
  stat_ecdf(geom = 'step', pad = F, size = 1) +
  scale_color_manual(values = cols) +
  labs(color = 'Riverine\nHW Severity') +
  xlab(bquote('Gross Primary Production (g'~O[2]~ m^-2~d^-1*')')) +
  ylab('ECDF') +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        legend.position = c(0.7,0.3),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

# width = 1500 height = 500
ggarrange(gpp_plot2,er_plot,nem_plot, ncol = 3)
