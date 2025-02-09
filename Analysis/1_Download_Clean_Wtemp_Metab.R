rm(list = ls())
dev.off()

library(tidyverse)
library(dataRetrieval)
library(daymetr)

# Metabolism Data Source
# https://www.sciencebase.gov/catalog/item/59bff507e4b091459a5e0982

setwd("D:/School/MichiganTech/Metabolism_Heatwave/Data/")

metab <- read_tsv("daily_predictions.tsv")
site <- read_tsv("site_data.tsv")

metab <- metab %>%
  mutate(GPP = round(GPP,2),
         ER = round(ER, 2)) %>%
  select(date, site_name, GPP, ER) %>%
  left_join(., site, by = 'site_name') %>%
  select(date, nwis_id, GPP, ER) %>%
  rename(site_no = nwis_id)

# Average number of metabolism estimates (in months) per site?
metab %>%
  group_by(site_no) %>%
  summarise(Count_GPP = sum(!is.na(GPP)),
            Count_ER = sum(!is.na(ER))) %>%
  ungroup() %>%
  summarise(Mean_GPP = mean(Count_GPP)/30,
            Median_GPP = median(Count_GPP)/30,
            SD_GPP = sd(Count_GPP)/30,
            Mean_ER = mean(Count_ER)/30,
            Median_ER = median(Count_ER)/30,
            SD_ER = sd(Count_ER)/30)

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
table(Wtemp_daily$Wtemp_cd)
Wtemp_daily <- Wtemp_daily %>%
  mutate(Wtemp = ifelse(Wtemp < 0, 0, Wtemp),
         Wtemp = ifelse(Wtemp_cd %in% c("A [4]","A <"), NA, Wtemp))

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
initial <- NROW(na.omit(Wtemp_daily$Wtemp)) # 317,614

Wtemp_daily <- Wtemp_daily %>%
  group_by(site_no) %>%
  arrange(-desc(Date)) %>%
  mutate(Wtemp_int = zoo::na.approx(Wtemp, maxgap = 2, na.rm=F),
         Wtemp = ifelse(is.na(Wtemp), Wtemp_int, Wtemp),) %>%
  select(!Wtemp_int) %>% 
  ungroup()

final <- NROW(na.omit(Wtemp_daily$Wtemp)) # 319,254

### How much water temp data was linearly interpolated? ----
round(((final - initial)/initial)*100, digits = 1) # 0.5%

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
  filter(!Date >= as.Date('2017-01-02')) %>%
  select(!c(tile, altitude, year, yday, site)) %>%
  relocate(site_no, Date, latitude, longitude)
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

initial2 <- NROW(na.omit(wmet$Wtemp)) # 295,869

wmet <- wmet %>%
  mutate(yhat = ifelse(yhat < 0,0,round(yhat,2)),
         corWtemp = ifelse(is.na(Wtemp), round(yhat,2), Wtemp))

final2 <- NROW(na.omit(wmet$corWtemp)) # 372,289

### How much water temp data was linearly interpolated? ----
round(((final2 - initial2)/initial2)*100, digits = 1) # 25.8%

wmet %>% 
  ggplot(aes(x = Wtemp, y = yhat)) +
  geom_point(alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = 'longdash', color = 'red') +
  stat_smooth(method = 'lm') +
  labs(x = 'Observed Water Temp (C)',
       y = 'Modeled Water Temp (C)') +
  theme_bw()

# Download daily mean discharge (Q) time series (need to remove tidal sites) ----
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
  distinct() %>%
  ungroup()

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

wtemp_discharge_metab <- left_join(wtemp_discharge, metab, by = c('site_no', 'date'))

# Linear interpolate discharge data gaps less than or equal to 2 days ----
initial3 <- NROW(na.omit(wtemp_discharge_metab$flow_cms)) # 291,411

wtemp_discharge_metab <- wtemp_discharge_metab %>%
  group_by(site_no) %>%
  arrange(-desc(date)) %>%
  mutate(
    flow_cms = zoo::na.approx(flow_cms, maxgap = 2, na.rm = FALSE)) %>%
  ungroup()

final3 <- NROW(na.omit(wtemp_discharge_metab$flow_cms)) # 291,411

wtemp_discharge_metab <- left_join(wtemp_discharge_metab, site, by = "site_no") %>%
  select(site_no,long_name,lat,lon,date,Atemp,Wtemp,flow_cms,GPP,ER)

write.csv(wtemp_discharge_metab,'wtemp_discharge_metab.csv', row.names = FALSE)

lat_long <- wtemp_discharge_metab %>%
  select(site_no, long_name, lat, lon) %>%
  distinct()

write.csv(lat_long, '48USGS_site_locations.csv', row.names = FALSE)
