##########################################################################
#######Effect of antibiotic duration on hospitalised patients ############
##################### Show withinhost changes    #########################
##########################################################################

# clean working environment
rm(list=ls()) 

# load libraries 
library(ggplot2)
library(ggpubr)
library(reshape)

# load data 
load('runs/populationgrowth_2021-10-26withinhost.Rdata')

# clean data 
low.s = dat[which(dat$s_growth < 0.25), c('day', 'S_short_perday', 'S_long_perday')]
low.s$growth_type = 'low.s'
colnames(low.s) = c('day', 'short_perday', 'long_perday', 'growth_type')
high.s = dat[which(dat$s_growth > 0.75), c('day', 'S_short_perday', 'S_long_perday')]
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
  dplyr::summarise_at(vars(short_perday, long_perday), mean)
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

dat.long$growth_type[grep('high', dat.long$growth_type)] = 'high'
dat.long$growth_type[grep('low', dat.long$growth_type)] = 'low'


dat.long$line.labs.x = 27
dat.long$line.labs = NA
dat.long$line.labs.y = NA

for (d in c(3, 15)){
  for (b in c('S', 'R')){
    for (g in c('high', 'low')){
      row = which(dat.long$abx.dur == d & dat.long$bact_type == b & dat.long$growth_type == g)
      dat.long$line.labs.y[row] = tail(dat.long$value[row], 1) + 0.025
      if (b == 'R' & g == 'low') {dat.long$line.labs.y[row] = tail(dat.long$value[row], 1) - 0.03}
      if (b == 'R' & g == 'low' & d == 3) {dat.long$line.labs.y[row] = tail(dat.long$value[row], 1) - 0.02}
      if (b == 'R' & g == 'high' & d ==15) {dat.long$line.labs.y[row] = tail(dat.long$value[row], 1) + 0.09}
      dat.long$line.labs[row] = paste0(ifelse(g == 'high', 'Fast', 'Slow'), 
                                       ' ', b, ' growth')
    }
  }
}

dat.long$arrow.right.seg.x = 0
dat.long$arrow.right.seg.x.end = 15
dat.long$arrow.right.seg.x.end[which(dat.long$dur_type == 'Short treatment duration')] = 3
dat.long$arrow.right.seg.y = 0.9
dat.long$arrow.right.seg.y.end = 0.9

dat.long$arrow.left.seg.x = 3
dat.long$arrow.left.seg.x.end = 0
dat.long$arrow.left.seg.y = 0.9
dat.long$arrow.left.seg.y.end = 0.9

dat.long$arrow.lab = 'Treatment duration'
dat.long$arrow.lab[which(dat.long$dur_type == 'Short treatment duration')] = 'Treatment\nduration'
dat.long$arrow.lab.y = 0.92
dat.long$arrow.lab.y[which(dat.long$dur_type == 'Short treatment duration')] = 0.94
dat.long$arrow.lab.x = 7.5
dat.long$arrow.lab.x [which(dat.long$dur_type == 'Short treatment duration')] = 1.5

ggplot(dat.long, aes(x = day, y = value, group = interaction(bact_type, growth_type), color = bact_type, linetype = growth_type)) + 
  geom_rect(data = dat.long[c(1, nrow(dat.long)),], aes(xmin = 0, xmax = abx.dur, ymin = -Inf, ymax = Inf, fill = 'grey', color = NA), alpha = 0.1) +
  geom_line(size = 1.5) + 
  #geom_point(alpha = 0.4) +
  #geom_smooth(size = 1.5, se = F, method = loess) +
  geom_text(aes(x = line.labs.x, y = line.labs.y, label = line.labs, color = bact_type), size = 5) +
  labs(y = 'Number of bacteria\n(Proportion of max gut carrying capacity)', x = 'Day since antibiotic started') +
  facet_grid(.~dur_type) + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 30)) +
  scale_color_manual(name = '', 
                     values = c('#D7816A', "#729B79"), 
                     labels = c('Resistant bacteria', 'Susceptible bacteria')) +
  scale_linetype(name = '', labels = c('Fast bacterial growth', 'Slow bacterial growth')) +
  geom_segment(aes(x = dat.long$arrow.right.seg.x, xend = dat.long$arrow.right.seg.x.end, 
                   y = dat.long$arrow.right.seg.y, yend = dat.long$arrow.right.seg.y.end),
               arrow = arrow(unit(4, "cm")),
               colour = "grey") +
  geom_segment(aes(x = dat.long$arrow.left.seg.x, xend = dat.long$arrow.left.seg.x.end, 
                   y = dat.long$arrow.left.seg.y, yend = dat.long$arrow.left.seg.y.end),
               arrow = arrow(unit(4, "cm")),
               colour = "grey") +
  geom_text(aes(x = arrow.lab.x, y = arrow.lab.y, label = arrow.lab), color = 'grey40', size = 4) +
  theme_minimal() +
  scale_fill_manual(values = c('grey20'), name = '', labels = 'Antibiotic treatment period') + 
  theme(legend.position = 'none', 
        strip.text = element_text(size = 15, face='bold'),
        text = element_text(size = 15), 
        legend.text = element_text(size=12), 
        legend.key = element_blank()) 

ggsave('~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/final_main/withinhost.png',  
       dpi = 500, 
       width = 13, height = 8)


