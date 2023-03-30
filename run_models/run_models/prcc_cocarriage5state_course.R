###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
############# explore PRCC of co-carriage 5-state model ###########################
###################################################################################

# NOTE: 
# this runs the mode_carriage5state_course function, this is different from the mode_carriage5state function in that:
# average number of days of antibiotics per ward per iteration is recorded 
# average number of courses of antibiotics per ward per iteration is recorded 
# duration of antibiotic treatment is drawn from a wider distribution (larger SD) so to compare the effect of duration vs number of courses 
# the outcome will be the number of R carriers per bed per day, and new acquisitions per admission (not the difference)

rm(list = ls()) # Clean working environment

################################### Dependencies and functions ################################################

# SAMPLE PARAMETER SPACE 
# load libraries 
require(pse) #load pse package for Latin Hypercube
require(sensitivity) #load sensitivity package for sensitivity analysis
require(parallel) # load parallel processing package to use multiple cores on computer (or cluster)

# record errors to debug 
cl <- makeCluster(detectCores(), outfile=paste0('error_files/parallel_error_cocarriage5state_course_', Sys.Date(), '.txt'))

# source functions on all cores
clusterCall(cl, function() {source('models/model_cocarriage5state_course.R')})
model = 'cocarriage5state_course'

# data.df is a dataframe of the parameter values in columns to feed into LHS function later
modelRun.cocarriage.course <- function (data.df) { 
  return(mapply(run_model_cocarriage5state_course, 
                data.df[,1], data.df[,2], data.df[,3], data.df[,4], data.df[,5], 
                data.df[,6], data.df[,7], data.df[,8], data.df[,9], 
                data.df[,10], data.df[,11], data.df[,12], data.df[,13], data.df[,14], data.df[,15], 
                data.df[,16], data.df[,17]
  ))
}

################################## Define parameters and run LHS ####################################
#list parameters, the probability density functions from which the parameter values will be calculated, and what are the arguments to these functions
#list parameters together with name, so they are "linked" or not easily confused
parameters <- list(
  n.bed = c("qunif", list(min=5, max=50)),            #n.bed; number of beds in the ward
  max.los = c("qunif", list(min=3, max=20)),          #max.los; mean of length of stay (exponential distribution)
  prop_R = c("qunif", list(min=0, max=0.8)),          #probability of initial carriage of resistant organisms
  prop_r = c("qunif", list(min=0, max=1)),            #proportion of S in (S+s): prob_start_S <- prop_S_nonR*(1-prob_R)
  prop_Sr = c("qunif", list(min = 0, max = 1)),           #proportion of Sr in (r+R): prob_start_Sr <- prop_Sr_inR*prob_R
  prop_S = c("qunif", list(min = 0, max = 1)),            #proportion of sr in (r+r): prob_start_sr <- prop_sr_inR*prob_R
  bif = c("qunif", list(min = 0, max = 1)),               #bacterial interference factor (pi_ssr = pi_r1 * bif )
  pi_ssr = c("qunif", list(min = 0.3, max = 0.6)),      #probability of being transmitted r to ss (ss—> ssr)
  repop.s = c("qunif", list(min = 0.02, max = 0.12)),     #probability of regrowth of S  (s—>S)
  repop.r = c("qunif", list(min = 0.02, max = 0.15)),     #probability of regrowth of s (sr—> sR)
  mu = c("qunif", list(min = 0.002, max=0.02)),         #probability of being decolonised to S (Sr—> S) 
  abx.s = c("qunif", list(min = 0.1, max=0.5)),         #probability of clearing S to become s
  abx.r = c("qunif", list(min = 0.1, max=0.5)),         #probability of clearing R to become r
  abx.type = c("qunif", list(min = 1, max = 1.0001)),
  p.infect = c("qunif", list(min = 0.1, max = 1)),        #probability of being prescribed narrow spectrum antibiotic
  cum.r.1 = c("qunif", list(min = 30, max=300)),        #admission day when cumulative probability of HAI requiring abx.r is 1
  meanDur = c("qunif", list(min = 3, max = 21))         #mean duration of antibiotics (normal distribution)
)

# arrange parameters in a way LHS will be happy with
q <- unlist(lapply(parameters, function(l) l[[1]]))
q.arg <- lapply(parameters, function(l) l[2:3])
factors <- names(parameters)

abxr.effectiveness = ifelse(parameters$abx.r$max < 0.01, 'zero', 'notzero')
N = 500

# Check order of variables - MAKE SURE the variable listing and ORDER MATCHES the variable listing input into run_model
source('models/model_cocarriage5state_course.R')
if(all(names(parameters) == formalArgs(run_model_cocarriage5state_course))){
  
  old = Sys.time() # get start time
  LHS.cocarriage = LHS(modelRun.cocarriage.course, factors, N = N, q, q.arg, nboot = 100, cl = cl) #N is the size of the hypercube
  print(Sys.time() - old) # print elapsed time difference in nice format
  image_name = paste0(model, "_", N, "_", abxr.effectiveness, "_", Sys.Date())
  save(LHS.cocarriage, file = paste0("runs/", image_name, ".Rdata"))
  
} else {
  
  stop("Error: Listing of parameters in `parameters` does not match that in run_model_cocarriage5state function.")
  
}

stopCluster(cl)

