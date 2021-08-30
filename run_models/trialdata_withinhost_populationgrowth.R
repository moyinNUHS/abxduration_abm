###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
############# get trial data of population growth model   #########################
###################################################################################

rm(list = ls()) # clean working environment

library(parallel)
cl <- makeCluster(detectCores(), outfile = paste0('error_files/parallel_error_populationgrowth_withinhost', Sys.time(), '.txt'))

################################### Dependencies and functions ################################################
source('models/get_output_trialdata_withinhost_populationgrowth.R')
clusterCall(cl, function() {source('models/get_output_trialdata_withinhost_populationgrowth.R')})

###################### Sample parameter values to run model #################################################

N = 100            # sets of parameters
iterations = 200    # per set of parameters
sampled.para = data.frame(
  n.bed = rep(40, N),                         # n.bed; in this instance refers to number of repeats per set of parameters 
  n.day = rep(31, N),                
  prop_R = runif(N, min = 0, max = 0.8),      # probability of initial carriage of resistant organisms
  r_thres = runif(N,min = 0.01, max = 0.2),   # r_thres = threshold amount of bacteria before R can be transmitted
  K = runif(N, min = exp(18), max = exp(24)), 
  total_prop = runif(N, min = 0.5, max = 0.9), # total proportion of holding capacity, the starting amount of enterobacteriaceae
  s_growth = runif(N, min = 0.001,max = 2),  # s_growth = amount transmitted on log scale
  fitness.r = runif(N, min = 0, max = 2), # fitness.r = growth constant for logistic growth
  r_trans = runif(N, min = 0.01, max = 0.5),  # r_trans = mean amount of R transmitted
  abx.s = runif(N, min = 0.1, max = 1),     # abx.s = amount of s killed by broad spectrum abx s
  abx.r = runif(N, min = 0.1, max = 1),         # abx.r = amount of r killed by broad spectrum abx r
  short_dur = rep(3, N),                      # mean short duration of narrow spectrum antibiotics 
  long_dur = rep(15, N),                      # mean long duration of narrow spectrum antibiotics 
  timestep = rep(1, N),                       
  iterations = rep(iterations, N)                         
)


########################################### Run the model #####################################################

if(all(names(sampled.para) == formalArgs(run_trialdata_withinhost_populationgrowth))){
  
  old = Sys.time() # get start time
  
  sampled.para.df.list = split(sampled.para, seq(nrow(sampled.para)))
  sampled.para.list = lapply(sampled.para.df.list, unlist)
  
  trial.data.populationgrowth = parLapply(cl, sampled.para.list, function(x) {
    
    print(x)
    
    trialdata = run_trialdata_withinhost_populationgrowth(n.bed = x['n.bed'], n.day = x['n.day'],
                                              prop_R  = x['prop_R'], r_thres = x['r_thres'],
                                              K = x["K"], total_prop = x["total_prop"], s_growth = x["s_growth"], fitness.r = x["fitness.r"], r_trans = x["r_trans"],
                                              abx.s  = x['abx.s'], abx.r  = x['abx.r'], 
                                              short_dur = x['short_dur'], long_dur = x['long_dur'], 
                                              timestep = x['timestep'], iterations = x['iterations'])
    trialdata$day = 0 : (nrow(trialdata) - 1)
    return(trialdata)
  })
  names(trial.data.populationgrowth) = paste0('sampled.para.row_', 1:nrow(sampled.para))
  print(Sys.time() - old) # print elapsed time difference in nice format
}

# combine into one dataset 
# get average of all iterations
dat.out = do.call('rbind', trial.data.populationgrowth)
dat = cbind(sampled.para[rep(1:nrow(sampled.para), each = nrow(trial.data.populationgrowth[[1]])),], dat.out)

# save data 
save(dat, file = paste0("runs/populationgrowth_", Sys.Date(), "_withinhost_botheff.Rdata"))

stopCluster(cl)

