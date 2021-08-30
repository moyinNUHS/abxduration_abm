###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
######### show proportion of R carriers with length of stay #######################
###################################################################################

rm(list = ls()) # clean working environment

library(ggplot2)
library(ggpubr)
library(grid)
library(dplyr)

# get data 
dat.list = lapply(paste0('runs/', list.files('runs')[grep('2021-05-19_varydur', list.files('runs'))]), function(x) get(load(x)))
dat = do.call('rbind', dat.list)

# clean data for plot 
### categorise data
dat$resistance.type = ifelse(dat$abx.r.abx.r < 0.01, 'cre', '3gcre')
dat$abx.used = 'Cephalosporin'
dat$abx.used[which(dat$p.r.day1.p.r.day1 > 0.5 & dat$p.r.after.p.r.after > 0.5)] = 'Carbapenem'
dat$abx.effect = 'Ineffective antibiotic'
dat$abx.effect[which(dat$resistance.type == '3gcre' & dat$abx.used == 'Carbapenem')] = 'Effective antibiotic'
### average outcomes per day 
plot.dat.wide = dat %>%
  group_by(days, abx.effect, dur.cat) %>%
  summarise_at(vars(newR, totalR), mean, na.rm = T)
plot.dat = reshape2::melt(plot.dat.wide, id.var = colnames(plot.dat.wide)[!colnames(plot.dat.wide) %in% c('newR', 'totalR')])
### order the factors so the correct order of graphs is produced
plot.dat$abx.effect = as.factor(plot.dat$abx.effect)
plot.dat$abx.effect = factor(plot.dat$abx.effect, 
                             levels = c('Ineffective antibiotic', 'Effective antibiotic'), 
                             labels = c('Ineffective antibiotic administered\nmajority of the time', 
                                        'Effective antibiotic administered\nmajority of the time'))
plot.dat$dur.cat = as.factor(plot.dat$dur.cat)
plot.dat$dur.cat = factor(plot.dat$dur.cat, 
                          levels = c('long_output', 'short_output'), 
                          labels = c('Long treatment duration', 
                                     'Short treatment duration'))
plot.dat$variable = as.factor(plot.dat$variable)
plot.dat$variable = factor(plot.dat$variable, 
                           levels = c('totalR', 'newR'), 
                           labels = c('Proportion of patients who are resistance\ncarriers on each hospitalisation day\n', 
                                      'Proportion of patients acquiring resistance carriage\nout of non-resistance carriers\non the previous hospitalisation day'))

plot.dat$baseline.lab = NA
plot.dat$baseline.lab[which(plot.dat$abx.effect == 'Effective antibiotic administered\nmajority of the time')] = "Average baseline prevalence of resistance\ncarriers admitted into the ward"

plot.dat$baseline.lab.x = NA
plot.dat$baseline.lab.x[which(plot.dat$abx.effect == 'Effective antibiotic administered\nmajority of the time')] = 15.5

plot.dat$baseline.lab.y = NA
plot.dat$baseline.lab.y[which(plot.dat$abx.effect == 'Effective antibiotic administered\nmajority of the time')] = 0.4


totalR = ggplot(data = plot.dat[which(plot.dat$variable == 'Proportion of patients who are resistance\ncarriers on each hospitalisation day\n'),], 
                aes(x = days, y = value, group = dur.cat, color = dur.cat)) + 
  geom_point(size = 1, alpha = 0.4) +
  geom_line() +
  #geom_smooth(method = 'lm', formula = y ~ splines::bs(x, 3), size = 1, alpha = 0.6, se = F) +
  geom_hline(yintercept = 0.4, linetype = 'dashed', color = 'red') + 
  labs(y = '', x = '') +
  geom_text(aes(x = baseline.lab.x, y = baseline.lab.y, label = baseline.lab), color = 'grey40', size = 4) +
  scale_color_manual(name = '', values = c('#D65780', '#93B5C6')) +
  scale_y_continuous(limits = c(0.1, 0.6), breaks = seq(0.1, 0.6, 0.1)) +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  theme_minimal() +
  facet_wrap(.~ abx.effect , ncol = 2, strip.position = "top") +
  theme(legend.position = 'none', 
        strip.text = element_text(size = 13),
        text = element_text(size = 13), 
        legend.text=element_text(size=13)) 
totalR.anno = annotate_figure(totalR, left = textGrob("Proportion of patients who were resistance\ncarriers on each hospitalisation day\n", 
                                   rot = 90, vjust = 0.5, gp = gpar(cex = 1.1)))

newR = ggplot(data = plot.dat[which(plot.dat$variable == 'Proportion of patients acquiring resistance carriage\nout of non-resistance carriers\non the previous hospitalisation day'),], 
                 aes(x = days, y = value, group = dur.cat, color = dur.cat)) + 
  geom_point(size = 1, alpha = 0.4) +
  geom_line() +
  #geom_smooth(method = 'lm', formula = y ~ splines::bs(x, 3), size = 1, alpha = 0.6, se = F) +
  labs(y = '', x = 'Day of hospitalisation') +
  scale_color_manual(name = '', values = c('#D65780', '#93B5C6')) +
  scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, 0.02)) +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  theme_minimal() +
  facet_wrap(. ~ abx.effect,  ncol = 2, strip.position = "left") +
  theme(legend.position = 'bottom', 
        strip.text.y = element_blank(),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 13), 
        legend.text=element_text(size=13))

newR.anno = annotate_figure(newR, left = textGrob("Proportion of patients who acquired resistance\ncarriage on each hospitalistion day out of\nnon-resistance carriers on the previous day", 
                                                       rot = 90, vjust = 0.5, gp = gpar(cex = 1.1)))

ggarrange(totalR.anno, newR.anno, nrow = 2, common.legend = T, legend = 'bottom')

ggsave('~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/final_main/Rwithdays.png',  
       dpi = 500,
       width = 9, height = 9)
