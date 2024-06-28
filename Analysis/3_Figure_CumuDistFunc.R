rm(list = ls())
dev.off()

library(tidyverse)
library(patchwork)

setwd("D:/School/MichiganTech/Metabolism_Heatwave/Data")

hw_metab <- read_csv('hw_metab.csv')

hw_metab <- hw_metab %>%
  mutate(abs_ER = abs(ER),
         category = fct_relevel(category, 'None', 'Moderate', 'Strong', 'Severe', 'Extreme'))

# Figures (ECDF) ----

cols = c("None" = 'blue',
         "Moderate" = '#FFC866' ,
         "Strong" = '#FF6900' ,
         "Severe" = '#9E0000',
         "Extreme" = '#2D0000')

(gpp_plot <- ggplot(data = hw_metab, aes(x = GPP, color = category)) +
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
          legend.text = element_text(size = 16)))

(er_plot <- ggplot(data = hw_metab, aes(x = abs_ER, color = category)) +
    stat_ecdf(geom = 'step', pad = F, linewidth = 1) +
    scale_color_manual(values = cols) +
    labs(x = expression(atop(Ecosystem~Respiration,(g~O[2]~m^-2~d^-1))),
         y = NULL) +
    scale_y_continuous(breaks = seq(0,1,0.2)) +
    coord_cartesian(xlim=c(0, 20)) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 16, color = "black"),
          axis.title.y = element_text(size = 16, color = "black"),
          axis.text.x = element_text(size = 16, color = "black"),
          axis.text.y = element_text(size = 16, color = "black"),
          legend.position = 'none'))

(nem_plot <- ggplot(data = hw_metab, aes(x = NEM, color = category)) +
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
          legend.position = 'none'))

# width = 1500 height = 500
gpp_plot + er_plot + nem_plot + plot_layout(nrow = 1)
