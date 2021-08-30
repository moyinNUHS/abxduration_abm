###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
############# get trial data of co-carriage 5-state model #########################
###################################################################################

# run the cocarriage model to vary DURATION 
# and get the patient matrix, abx matrix, carriage status matrix for each iteration

rm(list = ls()) # clean working environment

library(parallel)
cl <- makeCluster(detectCores(), outfile = paste0('error_files/parallel_error_cocarriage5state_varydur_', Sys.time(), '.txt'))

################################### Dependencies and functions ################################################
source('models/get_output_trialdata_varydur_cocarriage5state.R')
clusterCall(cl, function() {source('models/get_output_trialdata_varydur_cocarriage5state.R')})

###################### Sample parameter values (obtained from prcc) to run model ##############################
model.names = c('runs/cocarriage5state_400_zero_2021-05-17.Rdata', 
                'runs/cocarriage5state_400_notzero_2021-05-17.Rdata')

resistance.type = c('3gcre', 'cre')
abx.type = c('allbroadspectrum', 'allnarrowspectrum')

########################################### Run the model #####################################################
# run the model 4 times - cre vs 3gcre, all broad spectrum antibiotics vs all narrow spectrum antibiotics 

for (resist in resistance.type) {
  for (abx in abx.type){
    
    abxr.type = ifelse(resist == 'cre', '_zero', 'notzero')
    model.name = model.names[grep(abxr.type, model.names)]
    
    lhs.data = get(load(model.name))
    sampled.para = lhs.data$data
    
    sampled.para$max.los = runif(nrow(sampled.para), min=20, max=40)
    
    if (abx == 'allbroadspectrum') {
      sampled.para$p.r.day1 = 0.9 ; sampled.para$p.r.after = 0.9
    } else {
      sampled.para$p.r.day1 = 0.1 ; sampled.para$p.r.after = 0.1
    }
    
    if(all(names(sampled.para) == formalArgs(run_trialdata_varydur_cocarriage5state))){
      
      old = Sys.time() # get start time
      
      sampled.para.df.list = split(sampled.para, seq(nrow(sampled.para)))
      sampled.para.list = lapply(sampled.para.df.list, unlist)
      
      trial.data.cocarriage <- parLapply(cl, sampled.para.list, function(x) {
        run_trialdata_varydur_cocarriage5state(n.bed = x['n.bed'], max.los  = x['max.los'], 
                                               prop_R  = x['prop_R'], prop_r  = x['prop_r'], prop_Sr  = x['prop_Sr'], prop_S  = x['prop_S'],
                                               bif  = x['bif'], pi_ssr  = x['pi_ssr'], repop.s  = x['repop.s'], fitness.r  = x['fitness.r'],
                                               mu  = x['mu'], abx.s  = x['abx.s'], abx.r  = x['abx.r'], 
                                               p.infect  = x['p.infect'], cum.r.1  = x['cum.r.1'], p.r.day1  = x['p.r.day1'], p.r.after  = x['p.r.after'],
                                               short_dur  = x['short_dur'], long_dur  = x['long_dur'])
      })
      names(trial.data.cocarriage) = paste0('sampled.para.row_', 1:nrow(sampled.para))
      
      print(Sys.time() - old) # print elapsed time difference in nice format
    }
    
    # combine into one dataset 
    data.list = lapply(trial.data.cocarriage, function(one_paraset){
      
      days = rep(1 : length(one_paraset$short_output$prop_R_perday_perparameterset), 2)
      para.df = matrix(rep(one_paraset$para, each = length(days)), 
                       nrow = length(days))
      dur.cat = rep(names(one_paraset)[2:3], each = max(days))
      
      combine = data.frame(days = days, 
                           para = para.df, 
                           dur.cat = dur.cat, 
                           newR = c(one_paraset$short_output$prop_newR_perday_perparameterset, one_paraset$long_output$prop_newR_perday_perparameterset), 
                           totalR = c(one_paraset$short_output$prop_R_perday_perparameterset, one_paraset$long_output$prop_R_perday_perparameterset),
                           abxdur = c(rep(one_paraset$short_output$abxdur_perpatient_perparameterset, max(days)), 
                                      rep(one_paraset$long_output$abxdur_perpatient_perparameterset, max(days))), 
                           abx1dur = c(rep(one_paraset$short_output$abx1dur_perpatient_perparameterset, max(days)), 
                                       rep(one_paraset$long_output$abx1dur_perpatient_perparameterset, max(days))), 
                           abx2dur = c(rep(one_paraset$short_output$abx2dur_perpatient_perparameterset, max(days)), 
                                       rep(one_paraset$long_output$abx2dur_perpatient_perparameterset, max(days)))
      )
      colnames(combine)[grep('para', colnames(combine))] = names(one_paraset$para)
      
      return(combine)
      
    })
    dat = do.call('rbind', data.list)
    
    image_name = paste0(substr(model.name, 6, 26), resist, '_', abx, '_')
    save(dat, file = paste0("runs/", image_name, Sys.Date(), "_varydur.Rdata"))
    
  }
}

stopCluster(cl)






