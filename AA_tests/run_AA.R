#########################################################################
#######Effect of antibiotic duration on hospitalised patients############
################Determine number of iterations per run###################
#########################################################################

#=======================#
#       Set up          #
#=======================#

# clean working environment
rm(list = ls())

# set working directory 
setwd('~/Documents/nBox/git_projects/abxduration_abm')

# load libraries 
require(pse)         #load pse package for Latin Hypercube
require(sensitivity) #load sensitivity package for sensitivity analysis 
require(spartan)     #for AA 

#=======================#
#  Consistency testing  #
#=======================#
#resource: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002916

# NOTE: 
# run the models with iteration = 1 
iterationstotry = c(5, 50, 100, 125, 150)

torun = list(list(model = 'simple3state', abx.r.upper = 0.5, iterationstotry = iterationstotry), 
             list(model = 'simple3state', abx.r.upper = 0.00001, iterationstotry =  iterationstotry), 
             list(model = 'cocarriage5state', abx.r.upper = 0.5, iterationstotry =  iterationstotry), 
             list(model = 'cocarriage5state', abx.r.upper = 0.00001, iterationstotry =  iterationstotry), 
             list(model = 'populationgrowth', abx.r.upper = 0.5, iterationstotry =  iterationstotry), 
             list(model = 'populationgrowth', abx.r.upper = 0.00001, iterationstotry =  iterationstotry))

for (i in 1:length(torun)){
  
  dat = torun[[i]]
  model = dat[['model']]
  abx.r.upper = dat[['abx.r.upper']]
  iterationstotry = dat[['iterationstotry']]
  
  source('AA_tests/get_AA.R')
}

