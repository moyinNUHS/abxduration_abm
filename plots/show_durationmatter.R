###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
###################### show model exploration outputs #############################
###################################################################################

rm(list = ls()) # Clean working environment

#libraries
library(ggplot2)
library(dplyr)
library(purrr)
library(stringr)

#source plot function
source('plots/plot_functions/plot_durationmatter.R')

#load data
lhs.list.abs = list('runs/simple3state_420_notzero_2021-05-22.Rdata',
                    'runs/cocarriage5state_420_notzero_2021-05-23.Rdata',
                    'runs/populationgrowth_400_notzero_2021-05-24.Rdata',
                    'runs/simple3state_420_zero_2021-05-23.Rdata',
                    'runs/cocarriage5state_420_zero_2021-05-21.Rdata',
                    'runs/populationgrowth_400_zero_2021-05-24.Rdata')
lhs.list.treated = list('runs/simple3state_500_treated_notzero_2021-07-02.Rdata',
                        'runs/cocarriage5state_500_treated_notzero_2021-07-01.Rdata',
                        'runs/populationgrowth_500_treated_notzero_2021-07-02.Rdata',
                        'runs/simple3state_500_treated_zero_2021-07-01.Rdata',
                        'runs/cocarriage5state_500_treated_zero_2021-07-02.Rdata',
                        'runs/populationgrowth_500_treated_zero_2021-07-02.Rdata')

#put data into dataframe 
plot.data.totalR = make.plot.d(lapply(lhs.list.abs, make.df, outcome.type = 'totalR_diff'))
plot.data.newR = make.plot.d(lapply(lhs.list.abs, make.df, outcome.type = 'newR_diff'))
plot.data.treated = make.plot.d(lapply(lhs.list.treated, make.df, outcome.type = 'propR_treated_ward'))

plot.data = rbind.data.frame(plot.data.treated, plot.data.newR, plot.data.totalR)

#summarise data to get quantiles 
plot.d.list = list()
counter = 1
for (i in unique(as.character(plot.data$outcome.type))) {
  for (j in unique(as.character(plot.data$model.scenario))){
    
    d = plot.data[which(plot.data$model.scenario == j & plot.data$outcome.type == i),]
    
    plot.d.list[[counter]] = c(model = unique(as.character(d$model.name)),
                                        model.abxr = unique(as.character(d$model.abxr.type)), 
                                        scenario.outcome = paste(j, i),
                                        outcome = i,
                                        y = quantile(d$y, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = T))
    counter = counter + 1
  }
}
quant.data = do.call('rbind',  plot.d.list)
colnames(quant.data)[grep('y.', colnames(quant.data))] = c('2.5%', '25%', '50%', '75%', '97.5%')
quant.data = as.data.frame(quant.data)

quant.data$model = as.factor(quant.data$model)
quant.data$model = factor(quant.data$model, 
                          levels = c('Population growth', 'Co-carriage 5-state', 'Simple 3-state'))

treatedR_label = 'Resistance carriers amongst the treated patients'
totalR_label = 'Difference in proportion of resistance carriers per day\nbetween patients receiving long and short duration of antibiotics'
newR_label = 'Difference in proportion of non-resistance carriers acquiring resistance carriage during admission\nbetween patients receiving long and short duration of antibiotics'
quant.data$outcome = as.factor(quant.data$outcome)
quant.data$outcome = factor(quant.data$outcome, 
                            levels = c('propR_treated_ward','totalR_diff', 'newR_diff'),
                            labels = c(treatedR_label, totalR_label, newR_label))

quant.data$model.abxr = as.factor(quant.data$model.abxr)
quant.data$model.abxr = factor(quant.data$model.abxr, 
                               levels = c('notzero', 'zero'), 
                               labels = c('Treated with\neffective antibiotic', 
                                          'Not treated with\neffective antibiotic'))

quant.data$`2.5%` = as.numeric(quant.data$`2.5%`)
quant.data$`25%` = as.numeric(quant.data$`25%`)
quant.data$`50%` = as.numeric(quant.data$`50%`)
quant.data$`75%` = as.numeric(quant.data$`75%`)
quant.data$`97.5%` = as.numeric(quant.data$`97.5%`)

# add directional arrows 
quant.data$arrow.lab.left = NA
quant.data$arrow.lab.left[which(quant.data$outcome == totalR_label)] = 'More resistance carriers in the short duration ward\n'
quant.data$arrow.lab.right = NA
quant.data$arrow.lab.right[which(quant.data$outcome == totalR_label)] = 'More resistance carriers in the long duration ward\n'

arrow.lab.x = 0.25
arrow.lab.y = 2.4
arrow.y = 2.4
arrow.x = c(0, 0.5)

quant.data$arrow.lab.left.x = NA
quant.data$arrow.lab.left.x[which(quant.data$outcome == totalR_label)] = -arrow.lab.x 
quant.data$arrow.lab.right.x = NA
quant.data$arrow.lab.right.x[which(quant.data$outcome == totalR_label)] = arrow.lab.x 

quant.data$arrow.lab.left.y = NA
quant.data$arrow.lab.left.y[which(quant.data$outcome == totalR_label)] =  arrow.lab.y
quant.data$arrow.lab.right.y = NA
quant.data$arrow.lab.right.y[which(quant.data$outcome == totalR_label)] = arrow.lab.y

quant.data$arrow.lab.left.seg.x = NA
quant.data$arrow.lab.left.seg.x[which(quant.data$outcome == totalR_label)] = arrow.x[1]
quant.data$arrow.lab.left.seg.x.end = NA
quant.data$arrow.lab.left.seg.x.end[which(quant.data$outcome == totalR_label)] = -arrow.x[2]

quant.data$arrow.lab.right.seg.x = NA
quant.data$arrow.lab.right.seg.x[which(quant.data$outcome == totalR_label)] = arrow.x[1]
quant.data$arrow.lab.right.seg.x.end = NA
quant.data$arrow.lab.right.seg.x.end[which(quant.data$outcome == totalR_label)] = arrow.x[2]

quant.data$arrow.lab.left.seg.y = NA
quant.data$arrow.lab.left.seg.y[which(quant.data$outcome == totalR_label)] = arrow.y
quant.data$arrow.lab.left.seg.y.end = NA
quant.data$arrow.lab.left.seg.y.end[which(quant.data$outcome == totalR_label)] = arrow.y

quant.data$arrow.lab.right.seg.y = NA
quant.data$arrow.lab.right.seg.y[which(quant.data$outcome == totalR_label)] = arrow.y
quant.data$arrow.lab.right.seg.y.end = NA
quant.data$arrow.lab.right.seg.y.end[which(quant.data$outcome == totalR_label)] = arrow.y

quant.data$vline = NA 
quant.data$vline[which(quant.data$outcome != treatedR_label)] = 0

#plot graph 
# error bar plot to show distribution of the outputs
ggplot(quant.data, aes(y = model.abxr, color = model, group = model, x = `50%`)) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax =`97.5%`), position = position_dodge(width = 0.4), size = 2, height = 0) +
  geom_errorbarh(aes(xmin = `25%`, xmax =`75%`), position = position_dodge(width = 0.4), size = 4, height = 0) +
  geom_point(size = 3, colour = 'white', shape = 16, position = position_dodge(width = 0.4)) +
  facet_wrap(~ outcome , nrow = 3) +
  geom_vline(data = quant.data, aes(xintercept = vline), linetype="dotted", 
             color = "red", size = 1) +
  scale_color_manual(values = c('#15616D', '#70C1B3', '#FE938C'), name = 'Models') +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 16)) +
  scale_x_continuous(limit = (c(-arrow.x[2], max(quant.data$`97.5%`))), 
                     breaks = round(seq(round(-max(quant.data$`97.5%`) + 0.2, 1), round(max(quant.data$`97.5%`), 1), by = 0.1), 2)) + 
  labs(y = '', 
       x = 'Proportion of resistance carriers') +
  geom_text(aes(x = arrow.lab.left.x, y = arrow.lab.left.y,label = arrow.lab.left), color = 'black', size = 5) +
  geom_text(aes(x = arrow.lab.right.x, y = arrow.lab.right.y,label = arrow.lab.right), color = 'black', size = 5) +
  geom_segment(aes(x = arrow.lab.right.seg.x, xend = arrow.lab.right.seg.x.end, 
                   y = arrow.lab.right.seg.y, yend = arrow.lab.right.seg.y.end),
               arrow = arrow(unit(10, "cm")),
               colour = "grey") +
  geom_segment(aes(x = arrow.lab.left.seg.x, xend = arrow.lab.left.seg.x.end, 
                   y = arrow.lab.left.seg.y, yend = arrow.lab.left.seg.y.end),
               arrow = arrow(unit(10, "cm")),
               colour = "grey") +
  theme_bw() +
  theme(legend.position = 'bottom', 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = 'gray95'),
        text = element_text(size = 20)) + 
  guides(color = guide_legend(reverse = TRUE))


ggsave('~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/final_main/durationmatter.png',  
       width = 18, height = 14)

