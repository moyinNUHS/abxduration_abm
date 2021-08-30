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
working.directory.main = '~/Documents/nBox/git_projects/indiv_abxduration(inc metaanalysis)/indiv_abxduration/'
setwd(working.directory.main)

# load libraries 
require(pse)         #load pse package for Latin Hypercube
require(sensitivity) #load sensitivity package for sensitivity analysis 
require(spartan)     #for AA 

# workhorse function 
source('AA_tests/get_AA.R')

#=======================#
#  Consistency testing  #
#=======================#
#resource: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002916

# NOTE: 
# run the models with iteration = 1 
# make sure in the main working directory, AA_tests/AA_runs/`model`/abxrZERO and AA_tests/AA_runs/`model`/abxrNOTZERO already created

# simple 3 state model 
get_AA(model = 'simple3state', abx.r.upper = 0.5, iterationstotry =  c(5, 50, 100, 125, 150)) 
get_AA(model = 'simple3state', abx.r.upper = 0.00001, iterationstotry = c(5, 50, 100, 125, 150)) 

# cocarriage 5 state model 
get_AA(model = 'cocarriage5state', abx.r.upper = 0.5, iterationstotry = c(5, 50, 100, 125, 150)) 
get_AA(model = 'cocarriage5state', abx.r.upper = 0.00001, iterationstotry = c(5, 50, 100, 125, 150)) 

# population growth model 
get_AA(model = 'populationgrowth', abx.r.upper = 0.8, iterationstotry = c(5, 50, 100, 125, 150)) 
get_AA(model = 'populationgrowth', abx.r.upper = 0.00001, iterationstotry = c(5, 50, 100, 125, 150)) 

