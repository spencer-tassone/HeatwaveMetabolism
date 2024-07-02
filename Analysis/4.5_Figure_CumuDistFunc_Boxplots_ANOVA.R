rm(list = ls())
dev.off()

library(tidyverse)
library(egg)

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
    labs(color = 'Heatwave\nSeverity',
         x = expression(atop(Gross~Primary~Production,(g~O[2]~m^-2~d^-1))),
         y = 'Cumulative Distribution') +
    scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
    coord_cartesian(xlim = c(0, 15)) +
    scale_x_continuous(breaks = seq(0,15,5)) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20, color = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 20, color = "black"),
          legend.position = c(0.7,0.4),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20)))

(er_plot <- ggplot(data = hw_metab, aes(x = abs_ER, color = category)) +
    stat_ecdf(geom = 'step', pad = F, linewidth = 1) +
    scale_color_manual(values = cols) +
    labs(x = expression(atop(Ecosystem~Respiration,(g~O[2]~m^-2~d^-1))),
         y = NULL) +
    scale_y_continuous(breaks = seq(0,1,0.2)) +
    coord_cartesian(xlim=c(0, 20)) +
    scale_x_continuous(breaks = seq(0,20,5)) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = 'none'))

(nem_plot <- ggplot(data = hw_metab, aes(x = NEM, color = category)) +
    stat_ecdf(geom = 'step', pad = F, linewidth = 1) +
    scale_color_manual(values = cols) +
    labs(x = expression(atop(Net~Ecosystem~Metabolism,(g~O[2]~m^-2~d^-1))),
         y = NULL) +
    scale_y_continuous(breaks = seq(0,1,0.2)) +
    coord_cartesian(xlim=c(-16, 15)) +
    scale_x_continuous(breaks = seq(-15,15,5)) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = 'none'))


# HW Metabolism ANOVA ----
hw_metab %>%
  group_by(category) %>%
  summarise(Mean_GPP = mean(GPP, na.rm = TRUE),
            SD_GPP = sd(GPP, na.rm = TRUE),
            Mean_ER = mean(abs_ER, na.rm = TRUE),
            SD_ER = sd(abs_ER, na.rm = TRUE),
            Mean_NEM = mean(NEM, na.rm = TRUE),
            SD_NEM = sd(NEM, na.rm = TRUE))

aov_gpp <- aov(GPP~category, data = hw_metab)
aov_er <- aov(abs_ER~category, data = hw_metab)
aov_nem <- aov(NEM~category, data = hw_metab)

summary(aov_gpp) # p-value < 0.001
summary(aov_er) # p-value < 0.001
summary(aov_nem) # p-value < 0.001

tukey_gpp <- TukeyHSD(aov_gpp)
tukey_er <- TukeyHSD(aov_er)
tukey_nem <- TukeyHSD(aov_nem)

as.data.frame(tukey_gpp$category) %>%
  filter(`p adj` < 0.05) %>%
  mutate(`p adj` = round(`p adj`,5))

as.data.frame(tukey_er$category) %>%
  filter(`p adj` < 0.05) %>%
  mutate(`p adj` = round(`p adj`,5))

as.data.frame(tukey_nem$category) %>%
  filter(`p adj` < 0.05) %>%
  mutate(`p adj` = round(`p adj`,5))


# Figures (boxplot w/ sig) ----

(gpp_sig <- hw_metab %>%
   ggplot(aes(x = forcats::fct_rev(category), y = GPP, color = category)) +
   geom_boxplot(outlier.shape = NA) +
   stat_summary(fun = mean, geom = "point", shape = 4, size = 3) + # add the mean (line is median)
   coord_flip(ylim = c(0, 15)) +
   scale_y_continuous(breaks = seq(0,15,5)) +
   scale_color_manual(values = cols) +
   annotate('segment', x = c(5,4,4,3), xend = c(4,3,2,2), y = c(11,12,13,14), yend = c(11,12,13,14), color = 'black') +
   annotate('text', x = c(4.5,3.5,3,2.5), y = c(11.25,12.25,13.25,14.25),
            label = c('***', '***', '***', '*'),
            color = 'black', angle = -90, size = 8) +
   labs(y = expression(atop(Gross~Primary~Production,(g~O[2]~m^-2~d^-1))),
        x = 'Heatwave Severity') +
   theme_bw() +
   theme(axis.title.x = element_text(size = 20, color = "black"),
         axis.title.y = element_text(size = 20, color = "black"),
         axis.text.x = element_text(size = 20, color = "black"),
         axis.text.y = element_text(size = 20, color = "black"),
         legend.position = 'none',
         plot.margin = grid::unit(c(1,1,0,0), "mm")))

(er_sig <- hw_metab %>%
    ggplot(aes(x = forcats::fct_rev(category), y = abs_ER, color = category)) +
    geom_boxplot(outlier.shape = NA, aes(middle = mean(GPP))) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 3) + # add the mean (line is median)
    coord_flip(ylim = c(0,20)) +
    scale_y_continuous(breaks = seq(0,20,5)) +
    scale_color_manual(values = cols) +
    annotate('segment', x = c(5,5), xend = c(4,3), y = c(18,19), yend = c(18,19), color = 'black') +
    annotate('text', x = c(4.5,4), y = c(18.25,19.25),
             label = c('***','***'),
             color = 'black', angle = -90, size = 8) +
    labs(y = expression(atop(Ecosystem~Respiration,(g~O[2]~m^-2~d^-1))),
         x = NULL) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 20, color = "black"),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 20, color = "black"),
          axis.text.y = element_blank(),
          legend.position = 'none',
          plot.margin = grid::unit(c(1,1,0,0), "mm")))

(nem_sig <- hw_metab %>%
    ggplot(aes(x = forcats::fct_rev(category), y = NEM, color = category)) +
    geom_boxplot(outlier.shape = NA, aes(middle = mean(GPP))) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 3) + # add the mean (line is median)
    coord_flip(ylim = c(-16,15)) +
    scale_y_continuous(breaks = seq(-15,15,5)) +
    scale_color_manual(values = cols) +
    annotate('segment', x = c(5,5,4,4,3,3), xend = c(2,1,2,1,2,1), y = c(8,9.25,10.5,11.75,13,14.25), yend = c(8,9.25,10.5,11.75,13,14.25), color = 'black') +
    annotate('text', x = c(3.5,3,3,2.5,2.5,2), y = c(8.25,9.5,10.75,12,13.25,14.5),
             label = c('***','**','***', '**', '**', '**'),
             color = 'black', angle = -90, size = 8) +
    labs(y = expression(atop(Net~Ecosystem~Metabolism,(g~O[2]~m^-2~d^-1))),
         x = NULL) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 20, color = "black"),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 20, color = "black"),
          axis.text.y = element_blank(),
          legend.position = 'none',
          plot.margin = grid::unit(c(1,1,0,0), "mm")))

# width = 1500 height = 1000
ggarrange(gpp_plot,er_plot,nem_plot,gpp_sig,er_sig,nem_sig, nrow = 2)
