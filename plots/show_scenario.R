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
(seed = round(runif(1, 1, 10000))) #use 1280
seed = 1280
set.seed(seed)

# source code to run model
source("models/get_output_absdiff_cocarriage5state.R")

# course code to plot 
source('plots/plot_functions/plot_scenario.R')

#adjust default variables here
common.para = list(n.bed = 10,       # n.bed= number of beds in the ward
                   max.los = 7,     # mean.max.los= mean of max length of stay (exponential distribution)
                   n.day = 60,
                   short_dur = 5,
                   long_dur = 15,
                   meanDur = 7,
                   prop_R = 0.2,     # Probability of being colonized with resistant strain on admission
                   pi_ssr = 0.3,     # pi_ssr= probability of R transmitting 
                   p.infect.after = 0.1, 
                   p.infect = 0.1,   # p=probability of receiving antibiotic
                   p.r.day1 = 0.5, 
                   p.r.after = 0.5, 
                   timestep = 1)

cocarriage5state.para = list(prop_S = 0.5,      # Proportion of large S within non-resistant states (S+s)
                             prop_Sr = 0.4,                  
                             prop_r = 0.5,    
                             bif = 1,        # bacterial interference factor
                             mu = 0.2,        # mu= probability of clearance of Sr to become S
                             repop.s = 0.1, 
                             fitness.r = 1.1,    # probability of repopulation of Sr to become sR 
                             abx.s = 0.5,
                             abx.r = 0)
para = c(common.para, cocarriage5state.para)

dat = run_scenario(para, n.runs = 40)
  
plot_scenario(para, dat)

ggsave(paste0('~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/scenario.jpg'), 
       units = 'cm', width = 32, height=25)
