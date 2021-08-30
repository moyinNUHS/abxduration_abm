#########################################################################
#######Effect of antibiotic duration on hospitalised patients############
#############Show convergence of models over time########################
#########################################################################

# Clean working environment
rm(list = ls()) 

# load dependencies
source("models/los_abx_matrix_varydur.R")
source("models/msm_util_rtnorm.R")
source('plots/plot_functions/plot_stability.R')
source('plots/default_params.R')

# load libraries
library(reshape)
library(ggplot2)
library(ggpubr)

# iterations per run 
iterations = 50

# get data to plot 
simple_cre_data = get.stability.data(model = 'simple3state', iterations, scenario = 'cre')
simple_3gcre_data = get.stability.data(model = 'simple3state', iterations, scenario = '3gcre')

cocarriage5state_cre_data = get.stability.data(model = 'cocarriage5state', iterations, scenario = 'cre')
cocarriage5state_3gcre_data = get.stability.data(model = 'cocarriage5state', iterations, scenario = '3gcre')

populationgrowth_cre_data = get.stability.data(model = 'populationgrowth', iterations, scenario = 'cre')
populationgrowth_3gcre_data = get.stability.data(model = 'populationgrowth', iterations, scenario = '3gcre')

plot.data = rbind(simple_cre_data, simple_3gcre_data, 
                  cocarriage5state_cre_data, cocarriage5state_3gcre_data,
                  populationgrowth_cre_data, populationgrowth_3gcre_data)
plot.data$model = as.factor(plot.data$model)
plot.data$model = factor(plot.data$model, levels = c('simple3state', 'cocarriage5state', 'populationgrowth'), 
                         labels = c('Simple 3-state', 'Co-carriage 5-state', 'Population growth'))
plot.data$scenario = as.factor(plot.data$scenario)
plot.data$scenario = factor(plot.data$scenario, levels = c('3gcre', 'cre'), 
                            labels = c('Third-generation cephalosporin\nresistant Enterobacteriaceae',
                                       'Carbapenem-resistant Enterobacteriaceae'))

# save data 
save(plot.data, file = paste0('runs/stability_plot_data', Sys.Date(), '.Rdata'))

# plot 
# load('runs/stability_plot_data.Rdata')
p = ggplot(plot.data, aes(x = x, y = value, colour = as.factor(iter))) +
  geom_line(size = 0.4, alpha = 0.4) +
  scale_color_manual(values = rep('grey50', (iterations))) +
  labs(y = 'Difference in number of resistance carriers per day between\nthose receiving long vs short antibiotic treatment duration',
       x = 'Days of observation')+
  geom_vline(xintercept = 150, linetype='dashed', size = 0.7, colour = 'red')+
  theme_minimal() +
  facet_grid(model ~ scenario) + 
  theme(legend.position = 'none', 
        text = element_text(size = 15))

# save image as pdf 
png(file = "~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/stability.png", 
    res = 450, 
    width = 8, height = 8, units="in")
p
dev.off() 
