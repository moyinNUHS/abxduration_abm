###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
################# run simple 3-state model to get heatmap #########################
###################################################################################

rm(list = ls()) # Clean working environment

################################### Dependencies and functions ################################################

# load libraries 
require(parallel)    #load parallel processing package to use multiple cores on computer (or cluster)

############################################ Define parameters  ##############################################

N = 25           # sets of parameters
iterations = 10   # per set of parameters

prop_R = seq(0.01, 0.8, length.out = N)
pi_ssr = seq(0, 0.3, length.out = N)
samples = expand.grid(pi_ssr, prop_R)

############################################ Run ABS DIFF models  #############################################

# record errors to debug 
cl <- makeCluster(detectCores())

# source functions on all cores
clusterCall(cl, function() {source('models/get_output_absdiff_simple3state.R')})

for (abx.r.val in c(0, 0.5)){
    
    sampled.para = data.frame(
        n.bed = rep(30, N^2),                         # n.bed; in this instance refers to number of repeats per set of parameters 
        max.los = rep(10, N^2),                
        prop_R = samples$Var2,      # probability of initial carriage of resistant organisms
        prop_S = rep(0.5, N^2), 
        bif = rep(0.9, N^2),    
        pi_ssr = samples$Var1,  
        repop.s = rep(0.12, N^2),  # "repop.s" probability of ss repopulated to S (Palleja, Nature Biology, 2018 on gut recovery ~9 months)
        mu = rep(0.02, N^2),      # "mu", probability of decolonisation (Haggai Bar-Yoseph, JAC, 2016, decreasing colonization rates from 76.7% (95% CI=69.3%–82.8%) at 1 month to 35.2% (95% CI=28.2%–42.9%) at 12 months of follow-up)
        abx.s = rep(0.5, N^2),     # abx.s = amount of s killed by broad spectrum abx s
        abx.r = rep(abx.r.val, N^2),         # abx.r = amount of r killed by broad spectrum abx r
        p.infect = rep(0.5, N^2),     # "p.infect", probability of being prescribed antibiotics
        p.infect.after = rep(0.07, N^2),    # admission day when cummulative prabability of HAI requiring abx.r is 1
        p.r.day1 = rep(0.5, N^2),    # probability of being prescribed broad spectrum antibiotic on admission 
        p.r.after = rep(0.5, N^2),    # probability of being prescribed broad spectrum antibiotic after admission 
        short_dur = rep(3, N^2), 
        long_dur = rep(20, N^2),
        iterations = rep(iterations, N^2)                         
    )
    
    abxr.effectiveness = ifelse(unique(sampled.para$abx.r) < 0.01, 'zero', 'notzero')
    
    # Check order of variables - MAKE SURE the variable listing and ORDER MATCHES the variable listing input into run_model
    source('models/get_output_absdiff_simple3state.R')
    if(all(names(sampled.para)  == formalArgs(run_absdiff_simple3state))){
        
        old = Sys.time() # get start time
        
        sampled.para.df.list = split(sampled.para, seq(nrow(sampled.para)))
        sampled.para.list = lapply(sampled.para.df.list, unlist)
        
        out = parallel::parLapply(cl, sampled.para.list, function(x) {
            
            print(x)
            
            trialdata = run_absdiff_simple3state(n.bed= x['n.bed'], max.los = x['max.los'], 
                                                 prop_R= x['prop_R'], prop_S = x['prop_S'], 
                                                 bif = x['bif'], pi_ssr = x['pi_ssr'], repop.s = x['repop.s'], mu = x['mu'], 
                                                 abx.s = x['abx.s'], abx.r = x['abx.r'],
                                                 p.infect = x['p.infect'], p.infect.after = x['p.infect.after'], 
                                                 p.r.day1 = x['p.r.day1'], p.r.after = x['p.r.after'], 
                                                 short_dur = x['short_dur'], long_dur = x['long_dur'], 
                                                 iterations = x['iterations'])
            
            return(trialdata)
        })
        names(out) = paste0('sampled.para.row_', 1:nrow(sampled.para))
        print(Sys.time() - old) # print elapsed time difference in nice format
    }
    
    
    d = do.call('rbind', out)
    colnames(d) = res.names
    dat = cbind(sampled.para, d)
    
    save(dat, file = paste0("runs/simple3state_transheat_absdiff", abxr.effectiveness, '', Sys.Date(), ".Rdata"))
    
}

stopCluster(cl)

# ############################################ Run TREATED models  #############################################
# 
# # record errors to debug 
# cl <- makeCluster(detectCores(), outfile = paste0('error_files/parallel_error_simple3state_', Sys.Date(), '.txt'))
# 
# # source functions on all cores
# clusterCall(cl, function() {source('models/get_output_treated_simple3state.R')})
# 
# for (abx.r.val in c(0, 0.5)){
#     
#     sampled.para = data.frame(
#         n.bed = rep(30, N^2),                         # n.bed; in this instance refers to number of repeats per set of parameters 
#         max.los = rep(10, N^2),                
#         prop_R = samples$Var2,      # probability of initial carriage of resistant organisms
#         prop_S = rep(0.5, N^2), 
#         bif = rep(0.9, N^2),    
#         pi_ssr = samples$Var1,  
#         repop.s = rep(0.12, N^2),  # "repop.s" probability of ss repopulated to S (Palleja, Nature Biology, 2018 on gut recovery ~9 months)
#         mu = rep(0.02, N^2),      # "mu", probability of decolonisation (Haggai Bar-Yoseph, JAC, 2016, decreasing colonization rates from 76.7% (95% CI=69.3%–82.8%) at 1 month to 35.2% (95% CI=28.2%–42.9%) at 12 months of follow-up)
#         abx.s = rep(0.5, N^2),     # abx.s = amount of s killed by broad spectrum abx s
#         abx.r = rep(abx.r.val, N^2),         # abx.r = amount of r killed by broad spectrum abx r
#         p.infect = rep(0.5, N^2),     # "p.infect", probability of being prescribed antibiotics
#         cum.r.1 = rep(150, N^2),    # admission day when cummulative prabability of HAI requiring abx.r is 1
#         p.r.day1 = rep(0.5, N^2),    # probability of being prescribed broad spectrum antibiotic on admission 
#         p.r.after = rep(0.5, N^2),    # probability of being prescribed broad spectrum antibiotic after admission 
#         meanDur = rep(20, N^2),
#         iterations = rep(iterations, N^2)                         
#     )
#     
#     abxr.effectiveness = ifelse(unique(sampled.para$abx.r) < 0.01, 'zero', 'notzero')
#     
#     # Check order of variables - MAKE SURE the variable listing and ORDER MATCHES the variable listing input into run_model
#     source('models/get_output_treated_simple3state.R')
#     if(all(names(sampled.para)  == formalArgs(run_treated_simple3state))){
#         
#         old = Sys.time() # get start time
#         
#         sampled.para.df.list = split(sampled.para, seq(nrow(sampled.para)))
#         sampled.para.list = lapply(sampled.para.df.list, unlist)
#         
#         out = parallel::parLapply(cl, sampled.para.list, function(x) {
#             
#             print(x)
#             
#             trialdata = run_treated_simple3state(n.bed= x['n.bed'], max.los = x['max.los'], 
#                                                  prop_R= x['prop_R'], prop_S = x['prop_S'], 
#                                                  bif = x['bif'], pi_ssr = x['pi_ssr'], repop.s = x['repop.s'], mu = x['mu'], 
#                                                  abx.s = x['abx.s'], abx.r = x['abx.r'],
#                                                  p.infect = x['p.infect'], cum.r.1 = x['cum.r.1'], 
#                                                  p.r.day1 = x['p.r.day1'], p.r.after = x['p.r.after'], 
#                                                  meanDur = x['meanDur'],
#                                                  iterations = x['iterations'])
#             
#             return(trialdata)
#         })
#         names(out) = paste0('sampled.para.row_', 1:nrow(sampled.para))
#         print(Sys.time() - old) # print elapsed time difference in nice format
#     }
#     
#     
#     d = do.call('rbind', out)
#     colnames(d) = res.names
#     dat = cbind(sampled.para, d)
#     
#     save(dat, file = paste0("runs/simple3state_transheat_treated", abxr.effectiveness, '', Sys.Date(), ".Rdata"))
#     
# }
# 
# stopCluster(cl)
# 
