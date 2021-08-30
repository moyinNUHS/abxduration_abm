###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
############# get trial data of co-carriage 5-state model #########################
###################################################################################

# run the cocarriage model to vary COURSES
# and get the patient matrix, abx matrix, carriage status matrix for each iteration

rm(list = ls()) # clean working environment

library(parallel)
cl <- makeCluster(detectCores(), outfile = paste0('error_files/parallel_error_cocarriage5state_varystarttime_', Sys.time(), '.txt'))

################################### Dependencies and functions ################################################
source('models/get_output_trialdata_varystarttime_cocarriage5state.R')
clusterCall(cl, function() {source('models/get_output_trialdata_varystarttime_cocarriage5state.R')})

resistance.type = c('cre')
abx.type = c('allbroadspectrum', 'allnarrowspectrum')

###################### Sample parameter values to run model #################################################
# run the model 4 times - cre vs 3gcre, all broad spectrum antibiotics vs all narrow spectrum antibiotics 

for (resist in resistance.type) {
  for (abx in abx.type){
    
    abx.r.min = ifelse(resist == '3gcre', 0.1, 0)
    abx.r.max = ifelse(resist == '3gcre', 0.5, 0.00001)
    abx.type.min = ifelse(abx == 'allbroadspectrum', 2, 1)
    abx.type.max = ifelse(abx == 'allbroadspectrum', 2.000001, 1.00001)
    
    N = 400
    sampled.para = data.frame(
      n.bed = runif(N, min=5, max=50),            #n.bed; number of beds in the ward
      max.los = runif(N, min=3, max=20),          #max.los; mean of length of stay (exponential distribution)
      prop_R = runif(N, min=0, max=0.8),          #probability of initial carriage of resistant organisms
      prop_r = runif(N, min=0, max=1),            #proportion of S in (S+s): prob_start_S <- prop_S_nonR*(1-prob_R)
      prop_Sr = runif(N, min = 0, max = 1),           #proportion of Sr in (r+R): prob_start_Sr <- prop_Sr_inR*prob_R
      prop_S = runif(N, min = 0, max = 1),            #proportion of sr in (r+r): prob_start_sr <- prop_sr_inR*prob_R
      bif = runif(N, min = 0, max = 1),               #bacterial interference factor (pi_ssr = pi_r1 * bif )
      pi_ssr = runif(N, min = 0.3, max = 0.6),      #probability of being transmitted r to ss (ss—> ssr)
      repop.s = runif(N, min = 0.02, max = 0.12),     #probability of regrowth of S  (s—>S)
      fitness.r = runif(N, min = 0, max = 2),     #probability of regrowth of s (sr—> sR)
      mu = runif(N, min = 0.002, max=0.02),         #probability of being decolonised to S (Sr—> S) 
      abx.s = runif(N, min = 0.1, max=0.5),         #probability of clearing S to become s
      abx.r = runif(N, min = abx.r.min, max = abx.r.max),         #probability of clearing R to become r
      abx.type = runif(N, min = abx.type.min, max = abx.type.max),
      p.infect = runif(N, min = 0.1, max = 1),        #probability of being prescribed narrow spectrum antibiotic
      cum.r.1 = runif(N, min = 30, max=300),        #admission day when cumulative probability of HAI requiring abx.r is 1
      meanDur = runif(N, min = 3, max = 21)         #mean duration of antibiotics (normal distribution)
    )
    
    ########################################### Run the model #####################################################
    
    if(all(names(sampled.para) == formalArgs(run_trialdata_varystarttime_cocarriage5state))){
      
      old = Sys.time() # get start time
      
      sampled.para.df.list = split(sampled.para, seq(nrow(sampled.para)))
      sampled.para.list = lapply(sampled.para.df.list, unlist)
      
      trial.data.cocarriage = parLapply(cl, sampled.para.list, function(x) {
        run_trialdata_varystarttime_cocarriage5state(n.bed = x['n.bed'], max.los  = x['max.los'], 
                                                  prop_R  = x['prop_R'], prop_r  = x['prop_r'], prop_Sr  = x['prop_Sr'], prop_S  = x['prop_S'],
                                                  bif  = x['bif'], pi_ssr  = x['pi_ssr'], repop.s  = x['repop.s'], fitness.r  = x['fitness.r'],
                                                  mu  = x['mu'], abx.s  = x['abx.s'], abx.r  = x['abx.r'], abx.type = x['abx.type'],
                                                  p.infect  = x['p.infect'], cum.r.1  = x['cum.r.1'], 
                                                  meanDur  = x['meanDur'])
      })
      names(trial.data.cocarriage) = paste0('sampled.para.row_', 1:nrow(sampled.para))
      print(Sys.time() - old) # print elapsed time difference in nice format
    }
    
    # combine into one dataset 
    dat.out = do.call('rbind', trial.data.cocarriage)
    dat = cbind(sampled.para, dat.out)
    
    image_name = paste0('cocarriage', resist, '_', abx, '_')
    save(dat, file = paste0("runs/", image_name, Sys.Date(), "_varystarttime.Rdata"))
  }
}

stopCluster(cl)

