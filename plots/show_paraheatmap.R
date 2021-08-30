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
library(epiR)

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
                                   "p.infect", "p.r.day1", "cum.r.1", "p.r.after", #Number of antibiotic prescriptions 
                                   "pi_ssr", "r_thres", "r_trans","bif", #Transmission of resistant Enterobacteriaceae 
                                   "prop_S", "prop_Sr","prop_r","prop_R", "total_prop", "K", #Carriage status on day one of admission to the ward 
                                   "repop.s","s_growth", #Regrowth of susceptible Enterobacteriaceae 
                                   "fitness.r",  #Regrowth of resistant Enterobacteriaceae 
                                   "mu",#decolonisation 
                                   "abx.s", "abx.r", #antibiotics killing
                                   "short_dur", "long_dur")))
param.labels = c("n", "l",  #ward level
                 "omega_day1", "omega_day1.r", "omega_after",  "omega_after.r", #Number of antibiotic prescriptions 
                 "phi_s",  "Gamma", "tau", "b", #Transmission of resistant Enterobacteriaceae 
                 "p_S", "p_Sr","p_r","p_R", "rho_e", "K", #Carriage status on day one of admission to the ward 
                 "g_s","c_s", #Regrowth of susceptible Enterobacteriaceae 
                 "f",  #Regrowth of resistant Enterobacteriaceae 
                 "mu",#decolonisation 
                 "alpha_s", "alpha_r", #antibiotics killing
                 "t_short", "t_long")

#### Input run data required for plotting parameter heatmap 
lhs.list.abs = list('runs/simple3state_420_notzero_2021-05-22.Rdata',
                    'runs/cocarriage5state_420_notzero_2021-05-23.Rdata',
                    'runs/populationgrowth_400_notzero_2021-05-24.Rdata',
                    'runs/simple3state_420_zero_2021-05-23.Rdata',
                    'runs/cocarriage5state_420_zero_2021-05-21.Rdata',
                    'runs/populationgrowth_400_zero_2021-05-24.Rdata')
lhs.list.treated = list('runs/simple3state_500_treated_notzero_2021-07-02.Rdata',
                        'runs/cocarriage5state_500_treated_notzero_2021-07-01.Rdata',
                        'runs/populationgrowth_500_treated_notzero_2021-07-02.Rdata',
                        'runs/simple3state_500_treated_zero_2021-06-30.Rdata',
                        'runs/cocarriage5state_500_treated_zero_2021-06-30.Rdata',
                        'runs/populationgrowth_500_treated_zero_2021-07-02.Rdata')

### plot prcc for resistance carrier prevalence
# get ranking positions of the parameters 
ranks.list.totalR = lapply(lhs.list.abs, getposition, outcome.type = 'totalR_diff')
ranks.list.newR = lapply(lhs.list.abs, getposition, outcome.type = 'newR_diff')
ranks.list.treated = lapply(lhs.list.treated, getposition, outcome.type = 'propR_treated_delta_diffavg')

# prepare data for plot 
plot.df = clean_ranks_data(treated = ranks.list.treated, 
                           totalR = ranks.list.totalR, 
                           newR = ranks.list.newR)

png('~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/final_main/para_heatmap.png', 
       width = 9000, height = 5000, res = 500)
plot_paraheatmap(plot.df)
dev.off()
