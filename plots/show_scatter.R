###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
######## Plot scatter plot of different parameters to check monotonicity###########
###################################################################################

# clean working environment
rm(list = ls()) 

# load libraries 
library(ggplot2)
library(ggpubr)

# source function 
source('plots/plot_functions/plot_scatter.R')

# get plots
p = scatter(get(load('runs/simple3state_500_notzero_2021-10-30.Rdata')))
png(file = '~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/simple3state_scatter.png', 
    res = 400,
    width = 7, height = 7, units = 'in')
p
dev.off()

p = scatter(get(load('runs/cocarriage5state_500_notzero_2021-10-31.Rdata')))
png(file = '~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/cocarriage5state_scatter.png', 
    res = 400,
    width = 7, height = 7, units = 'in')
p
dev.off()

p = scatter(get(load('runs/populationgrowth_500_notzero_2021-10-31.Rdata'))) 
png(file = '~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/populationgrowth_scatter.png', 
    res = 400,
    width = 7, height = 8, units = 'in')
p
dev.off()
