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
working.directory.main = '~/Documents/nBox/git_projects/abxduration_abm/'
setwd(working.directory.main)

#=======================#
#      Get data         #
#=======================#

plot.list = list()
for (model in c('simple3state', 'cocarriage5state', 'populationgrowth')){
  setwd(paste0(working.directory.main, "AA_tests/AA_runs/", model))
  combine.d = rbind(read.csv('abxrNOTZERO/AA_ATestMaxAndMedians.csv'), 
                    read.csv('abxrZERO/AA_ATestMaxAndMedians.csv'))
  d = combine.d[, c('samplesize', 'totalRtreated_diffMaxA', 'totalR_diffMaxA', 'newR_diffMaxA')] # using max instead of median to be conservative
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
                         labels = c('Exclusive colonisation',
                                    'Co-colonisation', 
                                    'Within-host growth'))
plot.data$abxr = as.factor(plot.data$abxr)
plot.data$abxr = factor(plot.data$abxr, 
                        levels = c('3gcre', 'cre'), 
                        labels = c('Administered antibiotics active against\nboth susceptible and resistant organisms', 
                                   'Administered antibiotics active against\nonly resistant organisms'))

plot = ggplot(plot.data, aes(x = samplesize, y = value, group = scenario, color = variable)) + 
  geom_line() + 
  ylab('Maximum A-Test scores for each sample size') + 
  xlab('Sample size') +
  scale_color_manual(values = c('#545E56', '#23B5D3', '#BAFF29'), 
                     labels = c('Treated individuals',
                                'Overall ward population',
                                'Non-resistance carriers acquiring\nresistance during admission'), 
                     name = 'Model outputs in terms of the difference in proportion of\nresistance carriers between the long and short wards in:') +
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
        text = element_text(size = 15)) +
  guides(color = guide_legend(title.position="top"))

## plot 
png("~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/AA.png", 
    res = 400,
    width = 9, height = 10, units = 'in')
plot
dev.off() 


