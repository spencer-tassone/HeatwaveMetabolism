rm(list = ls())
dev.off()

library(tidyverse)
library(egg)

setwd("D:/School/MichiganTech/Metabolism_Heatwave/Data")

hw_metab <- read_csv('hw_metab.csv')

hw_metab <- hw_metab %>%
  mutate(abs_ER = abs(ER),
         category = fct_relevel(category, 'None', 'Moderate', 'Strong', 'Severe', 'Extreme')) %>% 
  group_by(site_no) %>% 
  mutate(z_score_gpp = round((GPP - mean(GPP, na.rm = TRUE)) / sd(GPP, na.rm = TRUE), 2),
         z_score_er = round((abs_ER - mean(abs_ER, na.rm = TRUE)) / sd(abs_ER, na.rm = TRUE), 2),
         z_score_nep = round((NEP - mean(NEP, na.rm = TRUE)) / sd(NEP, na.rm = TRUE), 2)) %>% 
  ungroup()

# What proportion of the time is NEM net heterotrophic (NEM < 0) during the baseline no HW condition?
hw_metab %>%
  filter(category == 'None') %>%
  select(NEP) %>% 
  summarise(total = sum(!is.na(NEP)),
            total_netHeterotrophic = sum(NEP < 0, na.rm = TRUE),
            prop_netHeterotrophic = (total_netHeterotrophic/total)*100)

# Percent change in metabolism during heatwaves ----
hw_metab_summary <- hw_metab %>%
  group_by(category) %>%
  summarise(
    Mean_GPP = round(mean(GPP, na.rm = TRUE),2),
    SD_GPP = round(sd(GPP, na.rm = TRUE),2),
    Mean_ER = round(mean(abs_ER, na.rm = TRUE),2),
    SD_ER = round(sd(abs_ER, na.rm = TRUE),2),
    Mean_NEP = round(mean(NEP, na.rm = TRUE),2),
    SD_NEP = round(sd(NEP, na.rm = TRUE),2)
  )

hw_metab_summary

# Calculate percent increase/decrease compared to the 'None' category
percent_change <- hw_metab_summary %>%
  mutate(
    Percent_Change_GPP = round((Mean_GPP - Mean_GPP[category == "None"]) / Mean_GPP[category == "None"] * 100, 1),
    Percent_Change_ER = round((Mean_ER - Mean_ER[category == "None"]) / Mean_ER[category == "None"] * 100, 1),
    Percent_Change_NEP = round((Mean_NEP - Mean_NEP[category == "None"]) / Mean_NEP[category == "None"] * 100, 1)
  ) %>% 
  select(category, Percent_Change_GPP, Percent_Change_ER, Percent_Change_NEP)

percent_change

# Kruskal-Wallis test for metabolism among heatwave severity class ----
kruskal.test(GPP~category, data = hw_metab) # p-value < 0.001
kruskal.test(abs_ER~category, data = hw_metab) # p-value < 0.001
kruskal.test(NEP~category, data = hw_metab) # p-value < 0.001

# Post-Hoc: Dunn's test with the Benjamini-Hochberg correction ----
dunn.test::dunn.test(hw_metab$GPP, hw_metab$category, kw = FALSE, method = 'bh', altp = TRUE, rmc = TRUE)
dunn.test::dunn.test(hw_metab$abs_ER, hw_metab$category, kw = FALSE, method = 'bh', altp = TRUE, rmc = TRUE)
dunn.test::dunn.test(hw_metab$NEP, hw_metab$category, kw = FALSE, method = 'bh', altp = TRUE, rmc = TRUE)

# Figures (boxplots) ----
cols = c("None" = 'blue',
         "Moderate" = '#FFC866' ,
         "Strong" = '#FF6900' ,
         "Severe" = '#9E0000',
         "Extreme" = '#2D0000')

(gpp_raw <- hw_metab %>%
   ggplot(aes(x = category, y = GPP, color = category)) +
   geom_boxplot(outlier.shape = NA, linewidth = 1, fatten = 1) +
   stat_summary(fun = mean, geom = "point", shape = 4, size = 3) + # add the mean (line is median)
   coord_cartesian(ylim = c(0, 10.5)) +
   scale_y_continuous(breaks = seq(0,10,2)) +
   scale_color_manual(values = cols) +
   labs(y = expression(atop(Metabolism~Rate,(g~O[2]~m^-2~d^-1))),
        x = NULL,
        title = 'GPP') +
   theme_bw() +
   theme(axis.title.y = element_text(size = 14, color = "black"),
         axis.text.x = element_blank(),
         axis.text.y = element_text(size = 14, color = "black"),
         legend.position = 'none',
         plot.title = element_text(size = 14, color = "black", hjust = 0.5),
         plot.margin = grid::unit(c(1,1,0,0), "mm")))

(er_raw <- hw_metab %>%
    ggplot(aes(x = category, y = abs_ER, color = category)) +
    geom_boxplot(outlier.shape = NA, linewidth = 1, fatten = 1) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 3) + # add the mean (line is median)
    coord_cartesian(ylim = c(0,20)) +
    scale_y_continuous(breaks = seq(0,20,5)) +
    scale_color_manual(values = cols) +
    labs(y = NULL,
         x = NULL,
         title = 'ER') +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14, color = "black"),
          legend.position = 'none',
          plot.title = element_text(size = 14, color = "black", hjust = 0.5),
          plot.margin = grid::unit(c(1,1,0,0), "mm")))

(nep_raw <- hw_metab %>%
    ggplot(aes(x = category, y = NEP, color = category)) +
    geom_boxplot(outlier.shape = NA, linewidth = 1, fatten = 1) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 3) + # add the mean (line is median)
    coord_cartesian(ylim = c(-16,10)) +
    scale_y_continuous(breaks = seq(-15,10,5)) +
    scale_color_manual(values = cols) +
    labs(y = NULL,
         x = NULL,
         title = 'NEP') +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14, color = "black"),
          legend.position = 'none',
          plot.title = element_text(size = 14, color = "black", hjust = 0.5),
          plot.margin = grid::unit(c(1,1,0,0), "mm")))

(gpp_zscore <- hw_metab %>%
    ggplot(aes(x = category, y = z_score_gpp, color = category)) +
    geom_boxplot(outlier.shape = NA, linewidth = 1, fatten = 1) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 3) + # add the mean (line is median)
    coord_cartesian(ylim = c(-3.25, 3.25)) +
    scale_y_continuous(breaks = seq(-3,3,1)) +
    scale_color_manual(values = cols) +
    labs(y = 'Z-score',
         x = NULL) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 14, color = "black"),
          axis.title.y = element_text(size = 14, color = "black"),
          axis.text.x = element_text(size = 14, color = "black"),
          axis.text.y = element_text(size = 14, color = "black"),
          legend.position = 'none',
          plot.margin = grid::unit(c(1,1,0,0), "mm")))

(er_zscore <- hw_metab %>%
    ggplot(aes(x = category, y = z_score_er, color = category)) +
    geom_boxplot(outlier.shape = NA, linewidth = 1, fatten = 1) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 3) + # add the mean (line is median)
    coord_cartesian(ylim = c(-3.25, 3.25)) +
    scale_y_continuous(breaks = seq(-3,3,1)) +
    scale_color_manual(values = cols) +
    labs(y = NULL,
         x = 'Heatwave Severity') +
    theme_bw() +
    theme(axis.title.x = element_text(size = 14, color = "black"),
          axis.title.y = element_text(size = 14, color = "black"),
          axis.text.x = element_text(size = 14, color = "black"),
          axis.text.y = element_text(size = 14, color = "black"),
          legend.position = 'none',
          plot.margin = grid::unit(c(1,1,0,0), "mm")))

(nep_zscore <- hw_metab %>%
    ggplot(aes(x = category, y = z_score_nep, color = category)) +
    geom_boxplot(outlier.shape = NA, linewidth = 1, fatten = 1) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 3) + # add the mean (line is median)
    coord_cartesian(ylim = c(-3.25, 3.25)) +
    scale_y_continuous(breaks = seq(-3,3,1)) +
    scale_color_manual(values = cols) +
    labs(y = NULL,
         x = NULL) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 14, color = "black"),
          axis.title.y = element_text(size = 14, color = "black"),
          axis.text.x = element_text(size = 14, color = "black"),
          axis.text.y = element_text(size = 14, color = "black"),
          legend.position = 'none',
          plot.margin = grid::unit(c(1,1,0,0), "mm")))

# width = 1300, height = 700
ggarrange(gpp_raw,er_raw,nep_raw,gpp_zscore,er_zscore,nep_zscore, nrow = 2)  
