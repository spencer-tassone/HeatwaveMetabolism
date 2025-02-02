rm(list = ls())
dev.off()

library(tidyverse)
library(truncnorm)
library(egg)

# Conceptual ECDF ----
# No heatwaves
norm <- data.frame(
  GPP = rtruncnorm(10000, a = 0, b = 15, mean = 2.5, sd = 3),
  ER = rtruncnorm(10000, a = 0, b = 20, 5, 5)) %>%
  mutate(NEP = GPP - abs(ER),
         category = 'None')

# Moderate strength heatwaves
mod <- data.frame(
  GPP = rtruncnorm(10000, a = 0, b = 15, mean = 4, sd = 3),
  ER = rtruncnorm(10000, a = 0, b = 20, 6.5, 5)) %>%
  mutate(NEP = GPP - abs(ER),
         category = 'Moderate')

# Strong strength heatwaves
strong <- data.frame(
  GPP = rtruncnorm(10000, a = 0, b = 15, mean = 2.5, sd = 3.1),
  ER = rtruncnorm(10000, a = 0, b = 20, 8, 5)) %>%
  mutate(NEP = GPP - abs(ER),
         category = 'Strong')

# Severe strength heatwaves
severe <- data.frame(
  GPP = rtruncnorm(10000, a = 0, b = 15, mean = 1.5, sd = 2.5),
  ER = rtruncnorm(10000, a = 0, b = 20, 9.5, 5)) %>%
  mutate(NEP = GPP - abs(ER),
         category = 'Severe')

# Extreme strength heatwaves
extreme <- data.frame(
  GPP = rtruncnorm(10000, a = 0, b = 15, mean = 0.5, sd = 1.75),
  ER = rtruncnorm(10000, a = 0, b = 20, 11.5, 5)) %>%
  mutate(NEP = GPP - abs(ER),
         category = 'Extreme')

# Creating conceptual figures:
conceptual_dat <- rbind(norm, mod, strong, severe, extreme)

cols = c("None" = 'blue',
         "Moderate" = '#FFC866' ,
         "Strong" = '#FF6900' ,
         "Severe" = '#9E0000',
         "Extreme" = '#2D0000')

conceptual_dat$category <- factor(conceptual_dat$category,
                            levels = c('None',
                                       'Moderate',
                                       'Strong',
                                       'Severe',
                                       'Extreme'))

(gpp_conceptual <- ggplot(data = conceptual_dat, aes(x = GPP, color = category)) +
    stat_ecdf(geom = 'step', pad = F, linewidth = 1.5) +
    scale_color_manual(values = cols) +
    labs(color = 'Riverine\nHW Severity',
         x = NULL,
         y = expression(atop(Cumulative~Distribution,(Conceptual)))) +
    scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
    coord_cartesian(xlim=c(0, 15)) +
    annotate('text', label = 'a)', x = 0.8, y = 0.95, size = 7) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, color = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14, color = "black"),
          legend.position = c(0.7,0.4),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14)))

(er_conceptual <- ggplot(data = conceptual_dat, aes(x = ER, color = category)) +
    stat_ecdf(geom = 'step', pad = F, linewidth = 1.5) +
    scale_color_manual(values = cols) +
    labs(x = NULL,
         y = NULL) +
    scale_y_continuous(breaks = seq(0,1,0.2)) +
    coord_cartesian(xlim=c(0, 20)) +
    annotate('text', label = 'b)', x = 1.2, y = 0.95, size = 7) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, color = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 4, color = "white"),
          legend.position = 'none'))

(nep_conceptual <- ggplot(data = conceptual_dat, aes(x = NEP, color = category)) +
    stat_ecdf(geom = 'step', pad = F, linewidth = 1.5) +
    scale_color_manual(values = cols) +
    labs(x = NULL,
         y = NULL) +
    coord_cartesian(xlim=c(-16, 17.5)) +
    scale_x_continuous(breaks = seq(-15,15,5)) +
    annotate('text', label = 'c)', x = -15, y = 0.95, size = 7) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, color = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 4, color = "white"),
          legend.position = 'none'))

# Actual ECDF ----
setwd("D:/School/MichiganTech/Metabolism_Heatwave/Data")

hw_metab <- read_csv('hw_metab.csv')

hw_metab <- hw_metab %>%
  mutate(abs_ER = abs(ER),
         category = fct_relevel(category, 'None', 'Moderate', 'Strong', 'Severe', 'Extreme'))

(gpp_actual <- ggplot(data = hw_metab, aes(x = GPP, color = category)) +
    stat_ecdf(geom = 'step', pad = F, linewidth = 1.5) +
    scale_color_manual(values = cols) +
    labs(color = 'Heatwave\nSeverity',
         x = expression(atop(Gross~Primary~Production,(g~O[2]~m^-2~d^-1))),
         y = expression(atop(Cumulative~Distribution,(Observed)))) +
    scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
    coord_cartesian(xlim = c(0, 15)) +
    scale_x_continuous(breaks = seq(0,15,5)) +
    annotate('text', label = 'd)', x = 2, y = 0.95, size = 7) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 14, color = "black"),
          axis.title.y = element_text(size = 14, color = "black"),
          axis.text.x = element_text(size = 14, color = "black"),
          axis.text.y = element_text(size = 14, color = "black"),
          legend.position = 'none',
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14)))

(er_actual <- ggplot(data = hw_metab, aes(x = abs_ER, color = category)) +
    stat_ecdf(geom = 'step', pad = F, linewidth = 1.5) +
    scale_color_manual(values = cols) +
    labs(x = expression(atop(Ecosystem~Respiration,(g~O[2]~m^-2~d^-1))),
         y = NULL) +
    scale_y_continuous(breaks = seq(0,1,0.2)) +
    coord_cartesian(xlim=c(0, 20)) +
    scale_x_continuous(breaks = seq(0,20,5)) +
    annotate('text', label = 'e)', x = 1.2, y = 0.95, size = 7) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 14, color = "black"),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 14, color = "black"),
          axis.text.y = element_blank(),
          legend.position = 'none'))

(nep_actual <- ggplot(data = hw_metab, aes(x = NEP, color = category)) +
    stat_ecdf(geom = 'step', pad = F, linewidth = 1.5) +
    scale_color_manual(values = cols) +
    labs(x = expression(atop(Net~Ecosystem~Production,(g~O[2]~m^-2~d^-1))),
         y = NULL) +
    scale_y_continuous(breaks = seq(0,1,0.2)) +
    coord_cartesian(xlim=c(-16, 17.5)) +
    scale_x_continuous(breaks = seq(-15,15,5)) +
    annotate('text', label = 'f)', x = -15, y = 0.95, size = 7) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 14, color = "black"),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 14, color = "black"),
          axis.text.y = element_blank(),
          legend.position = 'none'))

# width = 1300, height = 700
ggarrange(gpp_conceptual,er_conceptual,nep_conceptual,gpp_actual,er_actual,nep_actual, nrow = 2)   

