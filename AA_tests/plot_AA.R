###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
############################### Plot AA graphs ####################################
###################################################################################

#=======================#
#       Set up          #
#=======================#

# clean working environment
rm(list = ls())

# load library
library(ggplot2)

# set working directory 
working.directory.main = '~/Documents/nBox/git_projects/indiv_abxduration(inc metaanalysis)/indiv_abxduration/'
setwd(working.directory.main)

#=======================#
#      Get data         #
#=======================#

plot.list = list()
for (model in c('simple3state', 'cocarriage5state', 'populationgrowth')){
  setwd(paste0(working.directory.main, "AA_tests/AA_runs/", model))
  combine.d = rbind(read.csv('abxrNOTZERO/AA_ATestMaxAndMedians.csv'), 
                    read.csv('abxrZERO/AA_ATestMaxAndMedians.csv'))
  d = combine.d[, c('samplesize', 'totalR_diffMaxA', 'newR_diffMaxA')] # using max instead of median to be conservative
  d$abxr = rep(c('3gcre', 'cre'), each = nrow(d)/2)
  d$model = model
  d.melt = reshape2::melt(d, id.vars = c('samplesize', 'abxr', 'model'))
  plot.list[[model]] = d.melt
}
plot.data = do.call('rbind', plot.list)
plot.data$scenario = paste(plot.data$model, plot.data$abxr, plot.data$variable)
plot.data$model = as.factor(plot.data$model) 
plot.data$model = factor(plot.data$model, levels = c('simple3state', 
                                                     'cocarriage5state', 
                                                     'populationgrowth'),
                         labels = c('Simple 3-state',
                                    'Co-carriage 5-state', 
                                    'Population growth'))
plot.data$abxr = as.factor(plot.data$abxr)
plot.data$abxr = factor(plot.data$abxr, 
                        levels = c('3gcre', 'cre'), 
                        labels = c('Third-generation cephalosporin resistant', 'Carbapenem resistant'))

plot = ggplot(plot.data, aes(x = samplesize, y = value, group = scenario, color = variable)) + 
  geom_line() + 
  ylab('Maximum A-Test scores for each sample size') + 
  xlab('Sample size') +
  scale_color_manual(values = c('#545E56', '#23B5D3'), 
                     labels = c('Proportion of resistance carriers per day',
                                'Proportion of non-resistance carriers acquiring\nresistance during admission')) +
  geom_hline(yintercept = 0.56, linetype="dashed", color = "red", size = 0.3) +
  geom_hline(yintercept = 0.66, linetype="dashed", color = "red", size = 0.3) +
  geom_hline(yintercept = 0.73, linetype="dashed", color = "red", size = 0.3) +
  annotate("text", x = 130, y = 0.57, label = "Small difference", size = 3) + 
  annotate("text", x = 130, y = 0.67, label = "Medium difference", size = 3) + 
  annotate("text", x = 130, y = 0.74, label = "Large difference", size = 3) + 
  scale_x_continuous(breaks = seq(0, 150, by = 25)) +
  facet_grid(model ~ abxr) +
  theme_minimal() + 
  theme(legend.position = 'bottom', 
        legend.title = element_blank(), 
        text = element_text(size = 15))

## plot 
pdf("~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/AA.pdf", 
    width = 8, height = 10) 
plot
dev.off() 


