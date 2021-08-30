##########################################################################
#######Effect of antibiotic duration on hospitalised patients ############
##################### Show withinhost changes    ##########################
##########################################################################

# clean working environment
rm(list=ls()) 

# load libraries 
library(ggplot2)
library(ggpubr)
library(reshape)

# load data 
load('runs/populationgrowth_2021-05-25_withinhost.Rdata')

# clean data 
low.s = dat[which(dat$s_growth < 0.5), c('day', 'S_short_perday', 'S_long_perday')]
low.s$growth_type = 'low.s'
colnames(low.s) = c('day', 'short_perday', 'long_perday', 'growth_type')
high.s = dat[which(dat$s_growth > 1.5), c('day', 'S_short_perday', 'S_long_perday')]
high.s$growth_type = 'high.s'
colnames(high.s) = c('day', 'short_perday', 'long_perday', 'growth_type')
low.r = dat[which(dat$fitness.r < 0.5), c('day', 'R_short_perday', 'R_long_perday')]
low.r$growth_type = 'low.r'
colnames(low.r) = c('day', 'short_perday', 'long_perday', 'growth_type')
high.r = dat[which(dat$fitness.r > 1.5), c('day', 'R_short_perday', 'R_long_perday')]
high.r$growth_type = 'high.r'
colnames(high.r) = c('day', 'short_perday', 'long_perday', 'growth_type')

dat.cat = rbind.data.frame(low.s, high.s, low.r, high.r)

r_thres = mean(c(dat$r_thres_short, dat$r_thres_long))

dat.cut.avgbyday = dat.cat %>%
  group_by(day, growth_type) %>%
  summarise_at(vars(short_perday, long_perday), mean)
dat.long = reshape2::melt(dat.cut.avgbyday, id.vars = c('day', 'growth_type'))
dat.long$dur_type = 'short'
dat.long$dur_type[grep('long', dat.long$variable)] = 'long'

dat.long$dur_type = as.factor(dat.long$dur_type)
dat.long$dur_type = factor(dat.long$dur_type, 
                           labels = c('Long treatment duration', 'Short treatment duration'))

dat.long$bact_type = 'S'
dat.long$bact_type[grep('.r', dat.long$growth_type)] = 'R'

dat.long$abx.dur = unique(dat$short_dur)
dat.long$abx.dur[which(dat.long$dur_type == 'Long treatment duration')] = unique(dat$long_dur)

dat.long$thres_text = 'Threshold of resistant bacteria population size\nfor the individual to be a resistance carrier'
dat.long$thres_text[which(dat.long$dur_type == 'Long treatment duration')] = NA
dat.long$thres_text.x = 22
dat.long$thres_text.x[which(dat.long$dur_type == 'Long treatment duration')] = NA
dat.long$thres_text.y = 0.2
dat.long$thres_text.y[which(dat.long$dur_type == 'Long treatment duration')] = NA

dat.long$growth_type[grep('high', dat.long$growth_type)] = 'high'
dat.long$growth_type[grep('low', dat.long$growth_type)] = 'low'

ggplot(dat.long, aes(x = day, y = value, group = interaction(bact_type, growth_type), color = bact_type, linetype = growth_type)) + 
  geom_rect(data = dat.long[c(1, nrow(dat.long)),], aes(xmin = 0, xmax = abx.dur, ymin = -Inf, ymax = Inf, fill = 'grey', color = NA), alpha = 0.1) +
  # geom_line(size = 1.5) + 
  # geom_point(alpha = 0.4) +
  geom_smooth(size = 1.5, se = F, method = loess) +
  geom_hline(yintercept = r_thres, linetype = 'dotted', color = 'blue') +
  geom_text(aes(x = thres_text.x, y = thres_text.y, label = thres_text), color = 'grey40', size = 3.7) +
  labs(y = 'Proportion of bacteria out of the total gut carrying capacity', x = 'Day since antibiotic started') +
  facet_grid(.~dur_type) + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 30)) +
  scale_color_manual(name = '', 
                     values = c('#D7816A', "#729B79"), 
                     labels = c('Resistant bacteria', 'Susceptible bacteria')) +
  scale_linetype(name = '', labels = c('Fast bacterial growth', 'Slow bacterial growth')) +
  theme_minimal() +
  scale_fill_manual(values = c('grey20'), name = '', labels = 'Antibiotic treatment period') + 
  guides(color = guide_legend(override.aes = list(color = c('#D7816A', "#729B79"),
                                                  fill  = c("white", "white"), 
                                                  shape = c(15, 15),
                                                  size = c(5, 5),
                                                  linetype = 0)), 
         linetype = guide_legend(override.aes = list(fill  = c("white", "white")))) +
  theme(legend.position = 'bottom', 
        strip.text = element_text(size = 15),
        text = element_text(size = 15), 
        legend.text = element_text(size=12), 
        legend.key = element_blank()) 

ggsave('~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/final_main/withinhost.png',  
       dpi = 500, 
       width = 12, height = 8)


