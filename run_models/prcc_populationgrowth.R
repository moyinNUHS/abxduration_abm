###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
###################### explore PRCC of frequency model ############################
###################################################################################

#NOTE: 
# Prior to running the model, check that 
# 1. number of iterations for each set of parameters required - AA test 
# 2. number of days of burn-in required - plot equilibrium 

# After running the model, check that 
# 1. size of latin hypercube is adequate with SMBA (included below)

rm(list = ls()) # Clean working environment

################################### Dependencies and functions ################################################

# SAMPLE PARAMETER SPACE 
# load libraries 
require(pse) #load pse package for Latin Hypercube
require(sensitivity) #load sensitivity package for sensitivity analysis
require(parallel) # load parallel processing package to use multiple cores on computer (or cluster)

# record errors to debug 
cl <- makeCluster(detectCores(), outfile=paste0('error_files/parallel_error_populationgrowth_', Sys.Date(), '.txt'))

# source functions on all cores
clusterCall(cl, function() {source('models/get_output_absdiff_populationgrowth.R')})
model = 'populationgrowth'

modelRun.populationgrowth <- function (data.df) { #data.df is a dataframe of the parameter values in columns 
  return(mapply(run_absdiff_populationgrowth, 
                data.df[,1], data.df[,2], data.df[,3], 
                data.df[,4], data.df[,5], data.df[,6], 
                data.df[,7], data.df[,8], data.df[,9], 
                data.df[,10], data.df[,11], data.df[,12], 
                data.df[,13], data.df[,14], data.df[,15], 
                data.df[,16], data.df[,17], data.df[,18]
  ))
}

################################## Define parameters and run LHS ####################################
#list parameters, the probability density functions from which the parameter values will be calculated, and what are the arguments to these functions
#list parameters together with name, so they are "linked" or not easily confused
parameters <- list(
  n.bed = c("qunif", list(min = 5, max = 50)),         #n.bed; number of beds in the ward
  max.los = c("qunif", list(min = 3, max = 20)),       #max.los; mean of length of stay (exponential distribution)
  p.infect = c("qunif", list(min = 0.1, max = 1)),     # probability of being prescribed narrow spectrum antibiotic
  cum.r.1 = c("qunif", list(min = 30, max = 300)),     # admission day when cumulative prabability of HAI requiring abx.r is 1
  p.r.day1 = c("qunif", list(min = 0.1, max = 1)),     # probability of being prescribed broad spectrum antibiotic on day 1 of admission 
  p.r.after = c("qunif", list(min = 0.1, max = 1)),     # probability of being prescribed broad spectrum antibiotic after admission 
  K = c("qunif", list(min = exp(18), max = exp(24))),            # gut holding capacity, on log scale, largest R number possible is exp(300) - typical colonic bacteria 10^14 number/mL content https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4991899/
  total_prop = c("qunif", list(min = 0.1, max = 0.9)), # total proportion of holding capacity, the starting amount of enterobacteriaceae
  prop_R = c("qunif", list(min = 0, max = 0.8)),       # probability of a patient coming into the ward carrying R
  pi_ssr = c("qunif", list(min = 0.001,max = 0.3)),    # pi_ssr = daily probability of transmitting resistant E coli
  r_trans = c("qunif", list(min = 0.01, max = 0.5)),        # r_trans = mean amount of R transmitted
  fitness.r = c("qunif", list(min = 0.5, max = 2)),   # fitness.r = growth for R
  r_thres = c("qunif", list(min = 0.01, max = 0.2)),       # r_thres = threshold amount of bacteria before R can be transmitted
  s_growth = c("qunif", list(min = 0.1,max = 2)),   # s_growth = amount transmitted on log scale
  abx.s = c("qunif", list(min = 0.1, max = 1)),       # abx.s = amount of s killed by broad spectrum abx s
  abx.r = c("qunif", list(min = 0, max = 0.000001)),       # abx.r = amount of r killed by broad spectrum abx r
  short_dur = c("qunif", list(min = 3, max = 7)),      # mean short duration of narrow spectrum antibiotics (normal distribution) 
  long_dur = c("qunif", list(min = 14, max = 21))      # mean long duration of narrow spectrum antibiotics (normal distribution)
)

# arrange parameters in a way LHS will be happy with
q <- unlist(lapply(parameters, function(l) l[[1]]))
q.arg <- lapply(parameters, function(l) l[2:3])
factors <- names(parameters)

abxr.effectiveness = ifelse(parameters$abx.r$max < 0.01, 'zero', 'notzero')
N = 400

# Check order of variables - MAKE SURE the variable listing and ORDER MATCHES the variable listing input into run_model
source('models/get_output_absdiff_populationgrowth.R')
if(all(names(parameters) == formalArgs(run_absdiff_populationgrowth))){
  
  print('calculating 1st run')
  ## run 1
  old = Sys.time() # get start time
  LHS.populationgrowth = LHS(modelRun.populationgrowth, factors, N = N, q, q.arg, nboot = 100, cl = cl) #N is the size of the hypercube
  print(Sys.time() - old) # print elapsed time difference in nice format
  image_name = paste0(model, "_", N, "_", abxr.effectiveness, "_", Sys.Date())
  save(LHS.populationgrowth, file = paste0("runs/", image_name, ".Rdata"))
  
 # print('calculating 2nd run')
  
  # run 2
  # N = N + 20
  # old = Sys.time() # get start time
  # LHS.populationgrowth2 = LHS(modelRun.populationgrowth, factors, N = N, q, q.arg, nboot = 100, cl = cl)
  # print(Sys.time() - old) # print elapsed time difference in nice format
  # image_name = paste0(model, "_", N, "_", abxr.effectiveness, "_", Sys.Date())
  # save(LHS.populationgrowth2, file = paste0("runs/", image_name, ".Rdata"))
  # 
} else {
  
  stop("Error: Listing of parameters in `parameters` does not match that in run_absdiff_populationgrowth function.")
}

stopCluster(cl)

# Check agreement between runs to decide if our sample size ie size of LHS is adequate 
# Symmetric Best Measure of Agreement (SBMA) between the PRCC coefficients of two runs with different sample sizes.
(SBMA <- sbma(LHS.populationgrowth, LHS.populationgrowth2))
# value of -1 indicates complete disagreement between the runs 
# value of 1 indicated complete agreement  (>0.7 acceptable)
# caveat: if none of the model parameters is monotonically correlated with the output, 
# the agreement between runs may stay as low as 0.2 even for very large hypercubes.


