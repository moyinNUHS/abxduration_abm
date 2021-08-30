##########################################################################
#######Effect of antibiotic duration on hospitalised patients############
#####################Show one run per scenario ##########################
#########################################################################

# clean working environment
rm(list=ls()) 

# load libraries 
library(ggplot2)
library(ggpubr)
library(reshape)

# set seed for simulation
(seed = round(runif(1, 1, 10000))) #use 6245
seed = 6245
set.seed(seed)

# source code to run model
model = 'cocarriage5state'
source("models/get_output_absdiff_cocarriage5state.R")

# course code to plot 
source('plots/plot_functions/plot_scenario.R')

#adjust default variables here
common.para = list(esbl = list(n.bed = 20,       # n.bed= number of beds in the ward
                               max.los = 10,     # mean.max.los= mean of max length of stay (exponential distribution)
                               n.day = 300,
                               short_dur = 5,
                               long_dur = 15,
                               meanDur = 7,
                               prop_R = 0.8,     # Probability of being colonized with resistant strain on admission
                               pi_ssr = 0.2,     # pi_ssr= probability of R transmitting 
                               #simple cum.r.1 = 600, 
                               cum.r.1 = 900, 
                               p.infect = 0.5,   # p=probability of receiving antibiotic
                               p.r.day1 = 0.4, 
                               timestep = 1), 
                   cpe = list(n.bed = 10,       # n.bed= number of beds in the ward
                              max.los = 5,     # mean.max.los= mean of max length of stay (exponential distribution)
                              n.day = 30,
                              short_dur = 5,
                              long_dur = 15,
                              meanDur = 7,
                              prop_R = 0.2,     # Probability of being colonized with resistant strain on admission
                              pi_ssr = 0.05,     # pi_ssr= probability of R transmitting 
                              # simple cum.r.1 = 600, 
                              cum.r.1 = 300, 
                              p.infect = 0.3,   # p=probability of receiving antibiotic
                              p.r.day1 = 0.5, 
                              p.r.after = 0.5, 
                              timestep = 1))

cocarriage5state.para = list(esbl = list(prop_S = 0.5,     # Proportion of large S within non-resistant states (S+s)
                                         prop_Sr = 0.4,                  
                                         prop_r = 0.5,    
                                         bif = 0.7,        # bacterial interference factor
                                         mu = 0.01,        # mu= probability of clearance of Sr to become S
                                         repop.s = 0.01, 
                                         repop.r= 0.03,    # probability of repopulation of Sr to become sR 
                                         abx.s = 0.3,
                                         abx.r = 0.3),
                             cpe = list(prop_S = 0.5,      # Proportion of large S within non-resistant states (S+s)
                                        prop_Sr = 0.4,                  
                                        prop_r = 0.5,    
                                        bif = 0.7,        # bacterial interference factor
                                        mu = 0.2,        # mu= probability of clearance of Sr to become S
                                        repop.s = 0.1, 
                                        repop.r= 0.3,    # probability of repopulation of Sr to become sR 
                                        abx.s = 0.5,
                                        abx.r = 0))

plot_scenario(model = model, scenario = 'cpe')

ggsave(paste0('~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/scenario_', model,'.pdf'), 
       units = 'cm', width = 32, height=25)
