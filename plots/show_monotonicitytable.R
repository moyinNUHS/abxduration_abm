###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
#### Check monotonicity of the parameters - Hoeffding's D and Spearman's ##########
###################################################################################


# Clean working environment
rm(list = ls()) 

# Load libraries 
library(xtable)
library(Hmisc)
library(tidyverse)

# Source function to populate table 
source('plots/plot_functions/table_monotonicity.R')

# get tables 
raw.tab.simple = monotonicity.tab(get(load('runs/simple3state_420_notzero_2021-05-22.Rdata')))
#get.latex(raw.tab.simple, 'Simple 3-state model')


raw.tab.cocarriage = monotonicity.tab(get(load('runs/cocarriage5state_420_notzero_2021-05-23.Rdata')))
#get.latex(raw.tab.cocarriage, 'Cocarriage 5-state model')

raw.tab.populationgrowth = monotonicity.tab(get(load('runs/populationgrowth_400_notzero_2021-05-24.Rdata')))
#get.latex(raw.tab.populationgrowth, 'Population growth model')

combined = get.csv(raw.tab.simple, raw.tab.cocarriage, raw.tab.populationgrowth)
write.csv(combined, "runs/monotonicity_table.csv", row.names = FALSE)
