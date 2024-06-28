rm(list = ls())
dev.off()

library(tidyverse)
library(patchwork)
library(EflowStats)
library(lme4)
library(mgcv)

setwd("D:/School/MichiganTech/Metabolism_Heatwave/Data")

hw_metab <- read_csv('hw_metab.csv')

hw_metab <- hw_metab %>%
  mutate(abs_ER = abs(ER),
         category = fct_relevel(category, 'None', 'Moderate', 'Strong', 'Severe', 'Extreme'),
         site_no2 = as.factor(site_no),
         gpp_log = log10(GPP+1),
         er_log = log10(abs_ER+1),
         category2 = category,
         category2 = ifelse(category == 'None', 0, category2),
         category2 = ifelse(category == 'Moderate', 1, category2),
         category2 = ifelse(category == 'Strong', 2, category2),
         category2 = ifelse(category == 'Severe', 3, category2),
         category2 = ifelse(category == 'Extreme', 4, category2))

# Distributions
hist(hw_metab$i_max)
hist(hw_metab$GPP)
hist(hw_metab$gpp_log)
hist(hw_metab$abs_ER)
hist(hw_metab$er_log)

# Generalized additive models ----

hw_metab <- hw_metab %>%
  mutate(water_day = get_waterYearDay(date))

model_gpp <- gam(
  gpp_log ~ s(intensity_relThresh) + s(water_day, bs = 'cc') + s(site_no2, bs = 're'),
  data = hw_metab,
  family = gaussian
)

summary(model_gpp)
plot(model_gpp)

model_er <- gam(
  er_log ~ s(intensity_relThresh) + s(water_day, bs = 'cc') + s(site_no2, bs = 're'),
  data = hw_metab,
  family = gaussian
)

summary(model_er)
plot(model_er)

# Extract output of model 3 to remake figure in ggplot2
# Extract the smooth term components for s(intensity)
gpp_output <- plot.gam(model_gpp, pages = 1, seWithMean = TRUE)
er_output <- plot.gam(model_er, pages = 1, seWithMean = TRUE)

# Create a data frame with the values for intensity and the smooth term estimates
gpp_output_data <- data.frame(
  intensity = gpp_output[[1]]$x,
  intensity_fit = gpp_output[[1]]$fit,
  intensity_se.fit = gpp_output[[1]]$se,
  water_day = gpp_output[[2]]$x,
  water_day_fit = gpp_output[[2]]$fit,
  water_day_se.fit = gpp_output[[2]]$se,
  metab = 'GPP'
)

er_output_data <- data.frame(
  intensity = er_output[[1]]$x,
  intensity_fit = er_output[[1]]$fit,
  intensity_se.fit = er_output[[1]]$se,
  water_day = er_output[[2]]$x,
  water_day_fit = er_output[[2]]$fit,
  water_day_se.fit = er_output[[2]]$se,
  metab = 'ER'
)

# Add confidence intervals
gpp_output_data <- gpp_output_data %>%
  mutate(intensity_lower = intensity_fit - 2 * intensity_se.fit,
         intensity_upper = intensity_fit + 2 * intensity_se.fit,
         water_day_lower = water_day_fit - 2 * water_day_se.fit,
         water_day_upper = water_day_fit + 2 * water_day_se.fit)

er_output_data <- er_output_data %>%
  mutate(intensity_lower = intensity_fit - 2 * intensity_se.fit,
         intensity_upper = intensity_fit + 2 * intensity_se.fit,
         water_day_lower = water_day_fit - 2 * water_day_se.fit,
         water_day_upper = water_day_fit + 2 * water_day_se.fit)

output_data <- rbind(gpp_output_data, er_output_data)

# Function to convert water day of year value into water year date (2023 is a dummy year)
waterDOY_to_waterDATE <- function(day_of_year, start_year = 2023) {
  
  start_date <- as.Date(paste0(start_year, "-10-01"))
  
  date <- start_date + (day_of_year - 1)
  
  return(date)
}

output_data <- output_data %>%
  mutate(water_date = waterDOY_to_waterDATE(water_day))

(intensity_plot <- output_data %>%
    ggplot(aes(x = intensity, y = intensity_fit, group = metab)) +
    geom_line(aes(color = metab)) +
    geom_ribbon(aes(ymin = intensity_lower,
                    ymax = intensity_upper,
                    fill = metab),
                alpha = 0.4) +
    labs(x = expression(Heatwave~Intensity~(degree*C)),
         y = expression(Est.~Effect~of~s(Heatwave~Intensity))) +
    scale_x_continuous(breaks = seq(-2,8,2),
                       limits = c(-3.11,8.2)) +
    scale_y_continuous(breaks = seq(-0.5,0.5,0.1),
                       limits = c(-0.5,0.5),
                       labels = scales::label_number(accuracy = 0.1)) +
    scale_fill_manual(values = c("GPP" = "royalblue", "ER" = "orange"),
                      labels = c("GPP" = expression(log[10](GPP+1)),
                                 "ER" = expression(log[10](abs(ER)+1)))) +  
    scale_color_manual(values = c("GPP" = "royalblue", "ER" = "orange"),
                       labels = c("GPP" = expression(log[10](GPP+1)),
                                  "ER" = expression(log[10](abs(ER)+1)))) +
    annotate('text', x = 1, y = -0.2,
             label = expression(GPP~R[adj.]^2~"="~0.62),
             size = 5, hjust = 0) +
    annotate('text', x = 1, y = -0.35,
             label = expression(ER~R[adj.]^2~"="~0.53),
             size = 5, hjust = 0) +
    theme_bw() +
    theme(axis.text = element_text(color = 'black', size = 14),
          axis.title = element_text(color = 'black', size = 14),
          legend.position = c(0.2,0.8),
          legend.title = element_blank(),
          legend.text = element_text(color = 'black', size = 14)) +
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)))

(doy_plot <- output_data %>%
    ggplot(aes(x = water_date, y = water_day_fit, group = metab)) +
    geom_line(aes(color = metab)) +
    geom_ribbon(aes(ymin = water_day_lower,
                    ymax = water_day_upper,
                    fill = metab),
                alpha = 0.4) +
    labs(x = NULL,
         y = expression(Est.~Effect~of~s(Day~of~Year))) +
    scale_x_date(breaks = seq(as.Date('2023-10-01'),
                              as.Date('2024-09-30'),
                              '1 month'),
                 date_labels = "%b",
                 expand = c(0,0)) +
    scale_y_continuous(breaks = seq(-0.5,0.5,0.1),
                       limits = c(-0.5,0.5),
                       labels = scales::label_number(accuracy = 0.1)) +
    scale_fill_manual(values = c("GPP" = "royalblue", "ER" = "orange"),
                      labels = c("GPP" = expression(log[10](GPP+1)),
                                 "ER" = expression(log[10](abs(ER)+1)))) +  
    scale_color_manual(values = c("GPP" = "royalblue", "ER" = "orange"),
                       labels = c("GPP" = expression(log[10](GPP+1)),
                                  "ER" = expression(log[10](abs(ER)+1)))) +
    theme_bw() +
    theme(axis.text = element_text(color = 'black', size = 14),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title = element_text(color = 'black', size = 14),
          legend.position = 'none'))

# width = 600 height = 800
intensity_plot/doy_plot

# Linear mixed-effects models ----
# test_mod1 <- lmer(gpp_log ~ intensity_relThresh + (1 | site_no),
#                   data = hw_metab) # random intercept
# test_mod2 <- lmer(er_log ~ intensity_relThresh + (1 | site_no),
#                   data = hw_metab) # random intercept
# test_mod3 <- lmer(gpp_log ~ intensity_relThresh + (1 + intensity_relThresh | site_no),
#                   control = lmerControl(optimizer = "bobyqa",
#                                         optCtrl = list(maxfun = 100000)),
#                   data = hw_metab) # random slope and intercept
# test_mod4 <- lmer(er_log ~ intensity_relThresh + (1 + intensity_relThresh | site_no),
#                   control = lmerControl(optimizer = "bobyqa",
#                                         optCtrl = list(maxfun = 100000)),
#                   data = hw_metab) # random slope and intercept
# test_mod5 <- lmer(gpp_log ~ intensity_relThresh + category2 + (1 + intensity_relThresh | site_no),
#                   control = lmerControl(optimizer = "bobyqa",
#                                         optCtrl = list(maxfun = 100000)),
#                   data = hw_metab) # random slope and intercept
# test_mod6 <- lmer(er_log ~ intensity_relThresh + category2 + (1 + intensity_relThresh | site_no),
#                   control = lmerControl(optimizer = "bobyqa",
#                                         optCtrl = list(maxfun = 100000)),
#                   data = hw_metab) # random slope and intercept
# 
# summary(test_mod1)
# summary(test_mod2)
# summary(test_mod3)
# summary(test_mod4)
# summary(test_mod5)
# summary(test_mod6)
# 
# plot(test_mod1)
# qqnorm(resid(test_mod1))
# qqline(resid(test_mod1))
# 
# plot(test_mod2)
# qqnorm(resid(test_mod2))
# qqline(resid(test_mod2))
# 
# plot(test_mod3)
# qqnorm(resid(test_mod3))
# qqline(resid(test_mod3))
# 
# plot(test_mod4)
# qqnorm(resid(test_mod4))
# qqline(resid(test_mod4))
# 
# plot(test_mod5)
# qqnorm(resid(test_mod5))
# qqline(resid(test_mod5))
# 
# plot(test_mod6)
# qqnorm(resid(test_mod6))
# qqline(resid(test_mod6))
# 
# library(ggeffects)
# 
# pred.gpp <- ggpredict(test_mod1, terms = c('intensity'))
# pred.er <- ggpredict(test_mod2, terms = c('intensity'))
# 
# ggplot(pred.gpp) + 
#   geom_point(data = hw_metab,
#              aes(x = intensity, y = gpp_log, colour = category)) + 
#   geom_line(aes(x = x, y = predicted)) +          # slope
#   geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
#               fill = "lightgrey", alpha = 0.5) +  # error band
#   scale_color_manual(values = cols) +
#   theme_bw()
# 
# ggplot(pred.er) + 
#   geom_point(data = hw_metab,
#              aes(x = intensity, y = er_log, colour = category)) + 
#   geom_line(aes(x = x, y = predicted)) +          # slope
#   geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
#               fill = "lightgrey", alpha = 0.5) +  # error band
#   scale_color_manual(values = cols) +
#   theme_bw()
