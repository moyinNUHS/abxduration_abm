###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
####### Plot graphs to show the effect of no, low, high transmission on duration ##
###################################################################################

rm(list = ls()) # clean working environment

library(ggplot2)
library(ggpubr)
library(grid)
library(dplyr)

# load data 
# combine data and output of each iteration 
dat.list = lapply(paste0('runs/', list.files('runs/')[grep('2021-05-11_varytransmission', list.files('runs/'))]), function(x){
  
  lhs.data = get(load(x))
  lhs.data$abx.type = ifelse(length(grep('narrowspectrum', x)) > 0, 'Third-generation cephalosporin', 'Carbapenem')
  lhs.data$scenario = ifelse(all(lhs.data$abx.r.abx.r > 0.01), 'Third-generation cephalosporin resistance', 'Carbapenem resistance')
  return(lhs.data)
  
})
dat = do.call('rbind', dat.list)

plot.dat = dat[,c('days', 'sR', 'S', 's', 'sr', 'Sr', 'dur.cat', 'scenario', 'abx.type')]

# clean data for plot 

### order the factors so the correct order of graphs is produced
plot.dat$dur.cat = as.factor(plot.dat$dur.cat)
plot.dat$dur.cat = factor(plot.dat$dur.cat, 
                          levels = c('long_output', 'short_output'), 
                          labels = c('Long treatment duration', 
                                     'Short treatment duration'))

plot.dat.wide = plot.dat %>%
  group_by(days, dur.cat) %>%
  summarise_at(vars(sR, s, sr, Sr, S), mean)

plot.dat.long = reshape2::melt(plot.dat.wide, id.vars = c('days', 'dur.cat'))

plot.dat.long$variable = as.factor(plot.dat.long$variable)
plot.dat.long$variable = factor(plot.dat.long$variable, 
                                levels = c('sR', 'sr', 'Sr', 's', 'S'))

baseline = data.frame(dur.cat = "Short treatment duration", 
                      variable = NA, 
                      label = "Average baseline prevalence of resistance\ncarriers admitted into the ward", 
                      x = 300,
                      y = 0.375)
baseline$dur.cat = as.factor(baseline$dur.cat)
baseline$dur.cat = factor(baseline$dur.cat, 
                           levels = c('Long treatment duration', 
                                      'Short treatment duration'))

trans.rate = data.frame(dur.cat = rep(c("Long treatment duration", "Short treatment duration"), each = 3), 
                        variable = NA, 
                        label = rep(c("Transmission rate = 0", "Transmission rate = 0.1", "Transmission rate = 0.2"), 2), 
                        x = rep(c(75, 225, 375), 2),
                        y = rep(0.55, 6))

ggplot(plot.dat.long, aes(x = days, y = value, group = variable, color = variable)) + 
  geom_point(size = 0.5) + 
  geom_hline(yintercept = 0.4, linetype = 'dashed', color = 'red') +
  geom_line() + 
  geom_vline(xintercept = c(150, 300), linetype = 'dotted', color = 'grey') +
  labs(y = 'Proportion of patients', x = 'Day of observation in the wards') +
  geom_text(aes(x = x, y = y, label = label), color = 'gray40', size = 4, data = baseline) +
  geom_text(aes(x = x, y = y, label = label), color = 'black', size = 3, data = trans.rate) +
  scale_color_manual(name = 'Carriage type', values = c('#ee7a12', "#F2A359", "#f8caa0",  "#d6e1d9", "#AAC0AF")) +
  #scale_y_continuous(limits = c(0.2, 0.6), breaks = seq(0.2, 0.6, 0.1)) +
  scale_x_continuous(limits = c(0, 450), breaks = seq(0, 450, 50)) +
  theme_minimal() +
  facet_wrap(.~ dur.cat , ncol = 2, strip.position = "top") +
  theme(legend.position = 'bottom', 
        strip.text = element_text(size = 13),
        text = element_text(size = 13), 
        legend.text = element_text(size=13)) 

ggsave('~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/transmission_pergrp.png',  
       dpi = 500, 
       width = 10, height = 5)
