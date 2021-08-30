###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
################# explore PRCC of simple 3-state model ############################
###################################################################################

# NOTE: 
# Prior to running the model, check that 
# 1. number of iterations for each set of parameters required - AA test 
# 2. number of days of burn-in required - plot equilibrum 

# After running the model, check that 
# 1. size of latin hypercube is adequate with SMBA (included below)

rm(list = ls()) # Clean working environment

################################### Dependencies and functions ################################################

# load libraries 
require(pse)         #load pse package for Latin Hypercube
require(sensitivity) #load sensitivity package for sensitivity analysis 
require(parallel)    #load parallel processing package to use multiple cores on computer (or cluster)

# record errors to debug 
cl <- makeCluster(detectCores(), outfile = paste0('error_files/parallel_error_simple3state_treated_', Sys.Date(), '.txt'))

# source functions on all cores
clusterCall(cl, function() {source('models/get_output_treated_simple3state.R')})
model = 'simple3state'

# data.df is a dataframe of the parameter values in columns to feed into LHS function later
modelRun.simple <- function (data.df) { 
    return(mapply(run_treated_simple3state, 
                  data.df[,1], data.df[,2], data.df[,3], 
                  data.df[,4], data.df[,5], data.df[,6], 
                  data.df[,7], data.df[,8], data.df[,9], 
                  data.df[,10], data.df[,11], data.df[,12], 
                  data.df[,13],  data.df[,14], data.df[,15]
    ))
}

################################## Define parameters and run LHS ####################################
#list parameters, the probability density functions from which the parameter values will be calculated, and what are the arguments to these functions
#list parameters together with name, so they are "linked" or not easily confused
parameters <- list(
    n.bed = c("qunif", list(min=5, max=50)),         #"n.bed", number of beds in the ward
    max.los = c("qunif", list(min=3, max=20)),       #"max.los", mean of length of stay
    prop_R = c("qunif", list(min=0, max=0.8)),       #"prob_StartBact_R",probability of initial carriage of resistant organisms
    prop_S = c("qunif", list(min=0, max=1)),         #"prop_S", proportion of S in the population of S and ss
    bif = c("qunif", list(min=0, max=1)),            #"bif", bacterial interference factor
    pi_ssr = c("qunif", list(min=0.001, max=0.3)),   # "pi_ssr" probability of being transmitted r to ss (ss—> ssr)
    repop.s = c("qunif", list(min=0.02, max=0.12)),  # "repop.s" probability of ss repopulated to S (Palleja, Nature Biology, 2018 on gut recovery ~9 months)
    mu = c("qunif", list(min=0.002, max=0.02)),      # "mu", probability of decolonisation (Haggai Bar-Yoseph, JAC, 2016, decreasing colonization rates from 76.7% (95% CI=69.3%–82.8%) at 1 month to 35.2% (95% CI=28.2%–42.9%) at 12 months of follow-up)
    abx.s = c("qunif", list(min=0.1, max=0.5)),      # "abx.s", probability of S becoming ss after being on narrow spectrum antibiotics
    abx.r = c("qunif", list(min=0.1, max=0.5)),  # "abx.r", probability of R becoming ss after being on broad spectrum antibiotics
    p.infect = c("qunif", list(min=0.1, max=1)),     # "p.infect", probability of being prescribed antibiotics
    cum.r.1 = c("qunif", list(min=30, max= 300)),    # admission day when cummulative prabability of HAI requiring abx.r is 1
    p.r.day1 = c("qunif", list(min=0.1, max= 1)),    # probability of being prescribed broad spectrum antibiotic on admission 
    p.r.after = c("qunif", list(min=0.1, max= 1)),    # probability of being prescribed broad spectrum antibiotic after admission 
    meanDur = c("qunif", list(min=1, max=20))      # "meanDur", mean short duration of antibiotics (normal distribution)
)  

# arrange parameters in a way LHS will be happy with
q <- unlist(lapply(parameters, function(l) l[[1]]))
q.arg <- lapply(parameters, function(l) l[2:3])
factors <- names(parameters)

abxr.effectiveness = ifelse(parameters$abx.r$max < 0.01, 'zero', 'notzero')
N = 500

# Check order of variables - MAKE SURE the variable listing and ORDER MATCHES the variable listing input into run_model
source('models/get_output_treated_simple3state.R')
if(all(names(parameters) == formalArgs(run_treated_simple3state))){
    
    print('calculating 1st run')
    
    ## run 1 
    old = Sys.time() # get start time
    LHS.simple = LHS(modelRun.simple, factors, N = N, q, q.arg, nboot = 100, cl = cl) #N is the size of the hypercube
    print(Sys.time() - old) # print elapsed time difference in nice format
    image_name = paste0(model, "_", N, "_treated_", abxr.effectiveness, "_", Sys.Date())
    save(LHS.simple, file = paste0("runs/", image_name, ".Rdata"))
    
    # print('calculating 2nd run')
    # 
    # ## run 2 
    # N = N + 20
    # old = Sys.time() # get start time
    # LHS.simple2 = LHS(modelRun.simple, factors, N = N, q, q.arg, nboot = 100, cl = cl)
    # print(Sys.time() - old) # print elapsed time difference in nice format
    # image_name = paste0(model, "_", N, "_treated_", abxr.effectiveness, "_", Sys.Date())
    # save(LHS.simple2, file = paste0("runs/", image_name, ".Rdata"))
    
} else {
    
    stop("Error: Listing of parameters in `parameters` does not match that in run_treated_simple3state function.")
}

stopCluster(cl)

# Check agreement between runs to decide if our sample size ie size of LHS is adequate 
# Symmetric Best Measure of Agreement (SBMA) between the PRCC coefficients of two runs with different sample sizes.
(SBMA <- sbma(LHS.simple, LHS.simple2))
# value of -1 indicates complete disagreement between the runs 
# value of 1 indicated complete agreement  (>0.7 acceptable)
# caveat: if none of the model parameters is monotonically correlated with the output, 
# the agreement between runs may stay as low as 0.2 even for very large hypercubes.

