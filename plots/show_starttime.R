###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
# Plot graphs to show the effect of startlater courses vs one startonadmission course  #
###################################################################################

# clean working environment
rm(list=ls()) 

# load data 
# combine data and output of each iteration 
dat.list = lapply(paste0('runs/', list.files('runs/')[grep('25_varystarttime', list.files('runs/'))]), function(x){
  
  lhs.data = get(load(x))
  lhs.data$abx.type = ifelse(round(lhs.data$abx.type) == 1, 'Third-generation cephalosporin', 'Carbapenem')
  lhs.data$scenario = ifelse(all(lhs.data$abx.r > 0.1), 'Third-generation cephalosporin resistance', 'Carbapenem resistance')
  return(lhs.data)
  
})
dat = do.call('rbind', dat.list)


# clean data for plot 
admit.df = dat[, c('abxdur_startonadmission_perpatient_perparameterset', 
                   'prop_R_startonadmission_perparameterset',
                   'prop_newR_startonadmission_perparameterset',
                   'abx.type',
                   'scenario')]
colnames(admit.df) = c('dur', 'propR', 'newR','abx.type', 'scenario')
admit.df$intervention = 'startonadmission'
startlater.df = dat[, c('abxdur_startlater_perpatient_perparameterset', 
                        'prop_R_startlater_perparameterset',
                        'prop_newR_startlater_perparameterset',
                        'abx.type',
                        'scenario')]
colnames(startlater.df) = c('dur', 'propR', 'newR', 'abx.type', 'scenario')
startlater.df$intervention = 'startlater'

plot.dat.wide = rbind(admit.df, startlater.df)
rownames(plot.dat.wide) = NULL
plot.dat.wide = as.data.frame(apply(plot.dat.wide, 2, unlist))
plot.dat.wide$dur = as.numeric(plot.dat.wide$dur)
plot.dat.wide$propR = as.numeric(plot.dat.wide$propR)
plot.dat.wide$newR = as.numeric(plot.dat.wide$newR)
plot.dat = reshape2::melt(plot.dat.wide, id.var = colnames(plot.dat.wide)[!colnames(plot.dat.wide) %in% c('propR', 'newR')])

plot.dat$abx.effect = paste(plot.dat$scenario, plot.dat$abx.type)
plot.dat$abx.effect = ifelse(plot.dat$abx.effect == 'Third-generation cephalosporin resistance Carbapenem',
                             'Effective antibiotic administered\nmajority of the time', 
                             'Ineffective antibiotic administered\nmajority of the time')
plot.dat$abx.effect = factor(plot.dat$abx.effect, 
                             levels = c('Ineffective antibiotic administered\nmajority of the time', 
                                        'Effective antibiotic administered\nmajority of the time'))
plot.dat$scenario = as.factor(plot.dat$scenario)
plot.dat$scenario = factor(plot.dat$scenario, 
                           levels = c('Third-generation cephalosporin resistance', 
                                      'Carbapenem resistance'))
plot.dat$intervention = as.factor(plot.dat$intervention)
plot.dat$intervention = factor(plot.dat$intervention, 
                               levels = c('startonadmission', 'startlater'),
                               labels = c('Antibiotics when prescribed,\nwere started on admission', 
                                          'Antibiotics when prescribed,\nwere started later during admission'))
plot.dat$variable = as.factor(plot.dat$variable)
plot.dat$variable = factor(plot.dat$variable, 
                           levels = c('propR', 'newR'), 
                           labels = c('Proportion of patients who were\nresistance carriers per observation day', 
                                      'Proportion of patients who acquired\nresistance carriageduring admission'))

baseline = data.frame(variable = "Proportion of patients who were\nresistance carriers per observation day", 
                      abx.effect = "Effective antibiotic administered\nmajority of the time",
                      label = "Average baseline prevalence of resistance\ncarriers admitted into the ward", 
                      intervention = NA,
                      x = 12,
                      y = 0.4)
baseline$variable = as.factor(baseline$variable)
baseline$variable = factor(baseline$variable, 
                           levels = c('Proportion of patients who were\nresistance carriers per observation day', 
                                      'Proportion of patients who acquired\nresistance carriageduring admission'))
baseline$abx.effect = as.factor(baseline$abx.effect)
baseline$abx.effect = factor(baseline$abx.effect, 
                             levels = c('Effective antibiotic administered\nmajority of the time', 
                                        'Ineffective antibiotic administered\nmajority of the time'))


hline.y = data.frame(variable = c("Proportion of patients who were\nresistance carriers per observation day"), 
                     y = 0.4)
hline.y$variable = as.factor(hline.y$variable)
hline.y$variable = factor(hline.y$variable, 
                          levels = c('Proportion of patients who were\nresistance carriers per observation day', 
                                     'Proportion of patients who acquired\nresistance carriage during admission'))

# plot and save
p = ggplot(plot.dat, aes(x = dur, y = value, group = intervention, color = intervention)) +
  geom_point(alpha = 0.2, size = 0.2) + 
  geom_line(stat = "smooth", method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, alpha = 0.8, size = 1) +
  geom_text(data = baseline, aes(x = x, y = y, label = label), color = 'gray40', size = 4) +
  labs(x = 'Mean days of antibiotic treatment per patient', 
       y = '') +
  scale_x_continuous(limits = c(0, 20)) +
  scale_color_manual(values = c('#D62246', '#17BEBB'), name = '') + 
  facet_grid(variable ~ abx.effect, switch = "y") + 
  geom_hline(data = hline.y, aes(yintercept = y), linetype = 'dashed', color = 'red') + 
  theme_minimal() + 
  theme(legend.position = 'bottom', 
        text = element_text(size = 15))

ggsave(file = "~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/starttime.png", 
       dpi = 500, 
       width = 8.5, height = 8.5)

