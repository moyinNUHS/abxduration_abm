###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
######### Plot heatmap of the parameters against models ###########################
###################################################################################

rm(list = ls()) # Clean working environment

####Plot a heatmap of the parameters against models 
#libraries
library(reshape)
library(plyr)
library(scales)
library(ggplot2); library(ggpubr)
library(epiR); library(grid)

####Source models to get parameters to plot
source('models/get_output_absdiff_simple3state.R')
source('models/get_output_absdiff_cocarriage5state.R')
source('models/get_output_absdiff_populationgrowth.R')
source('plots/plot_functions/plot_paraheatmap.R')

# categorise the parameters
parameters = factor(unique(c(names(unlist(as.list(args(run_absdiff_simple3state)))), 
                             names(unlist(as.list(args(run_absdiff_cocarriage5state)))), 
                             names(unlist(as.list(args(run_absdiff_populationgrowth)))))), 
                    levels = rev(c("n.bed", "max.los",  #ward level
                                   "p.infect", "p.r.day1", "p.infect.after", "p.r.after", #Number of antibiotic prescriptions 
                                   "prop_S", "prop_Sr","prop_r","prop_R", "total_prop", "K", #Carriage status on day one of admission to the ward 
                                   "repop.s","s_growth", #Regrowth of susceptible Enterobacteriaceae 
                                   "fitness.r",  #Regrowth of resistant Enterobacteriaceae 
                                   "mu",#decolonisation 
                                   "abx.s", "abx.r", #antibiotics killing
                                   "pi_ssr", "r_thres", "r_trans","bif", #Transmission of resistant Enterobacteriaceae 
                                   "short_dur", "long_dur")))
param.labels = c("n", "l",  #ward level
                 "omega_day1", "omega_day1.r", "omega_after",  "omega_after.r", #Number of antibiotic prescriptions 
                 "p_S", "p_Sr","p_r","p_R", "rho_e", "K", #Carriage status on day one of admission to the ward 
                 "g_s",  "c_s", #Regrowth of susceptible Enterobacteriaceae 
                 "f",  #Regrowth of resistant Enterobacteriaceae 
                 "mu",#decolonisation 
                 "alpha_s", "alpha_r", #antibiotics killing
                 "phi_s",  "Gamma", "tau", "b", #Transmission of resistant Enterobacteriaceae 
                 "t_short", "t_long")

#### Input run data required for plotting parameter heatmap 
lhs.list = list('runs/simple3state_500_notzero_2021-10-30.Rdata',
                'runs/cocarriage5state_500_notzero_2021-10-31.Rdata',
                'runs/populationgrowth_500_notzero_2021-10-31.Rdata',
                'runs/simple3state_500_zero_2021-10-30.Rdata',
                'runs/cocarriage5state_500_zero_2021-10-31.Rdata',
                'runs/populationgrowth_500_zero_2021-10-31.Rdata')

### plot prcc for resistance carrier prevalence
# get ranking positions of the parameters 
ranks.list.totalR = lapply(lhs.list, getposition, outcome.type = 'totalR_diff')
ranks.list.newR = lapply(lhs.list, getposition, outcome.type = 'newR_diff')
ranks.list.treated = lapply(lhs.list, getposition, outcome.type = 'totalRtreated_diff')

# prepare data for plot 
plot.df = clean_ranks_data(treated = ranks.list.treated, 
                           totalR = ranks.list.totalR, 
                           newR = ranks.list.newR)

p = plot_paraheatmap(plot.df)

png('~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/final_main/para_heatmap.png', 
    width = 8500, height = 5000, res = 500)
plot_paraheatmap(plot.df)
dev.off()
