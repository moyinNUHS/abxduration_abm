###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
################# run co-carriage 5-state model to get heatmap #########################
###################################################################################

rm(list = ls()) # Clean working environment

################################### Dependencies and functions ################################################

# load libraries 
require(pse)         #load pse package for Latin Hypercube
require(parallel)    #load parallel processing package to use multiple cores on computer (or cluster)

############################################ Define parameters  ##############################################

N = 15           # sets of parameters
iterations = 10   # per set of parameters

meanDur = seq(0, 20, length.out = N)
prop_R = seq(0.01, 0.8, length.out = N)
pi_ssr = seq(0, 0.3, length.out = N)

samples.t = cbind(expand.grid(meanDur, pi_ssr), rep(0.1, N^2))
colnames(samples.t) = c('meanDur', 'pi_ssr', 'prop_R')
samples.p = cbind(expand.grid(meanDur, prop_R), rep(0.3, N^2))
colnames(samples.p) = c('meanDur', 'prop_R',  'pi_ssr')

############################################ Run ABS DIFF models  #############################################

# record errors to debug 
cl <- makeCluster(detectCores())

# source functions on all cores
clusterCall(cl, function() {source('models/get_output_abs_cocarriage5state.R')})

dat = list()

for (p in c('pi_ssr', 'prop_R')){ #parameter to vary
    
    if (p == 'pi_ssr'){
        samples = samples.t
    } else {
        samples = samples.p
    }
    
    sampled.para = data.frame(
        n.bed = rep(30, N^2),                         # n.bed; in this instance refers to number of repeats per set of parameters 
        max.los = rep(10, N^2),                
        prop_R = samples[['prop_R']],      # probability of initial carriage of resistant organisms
        prop_r = rep(0.5, N^2), 
        prop_Sr = rep(0.5, N^2),
        prop_S = rep(0.5, N^2), 
        bif = rep(0.9, N^2),    
        pi_ssr = samples[['pi_ssr']],  
        repop.s = rep(0.12, N^2),  # "repop.s" probability of ss repopulated to S (Palleja, Nature Biology, 2018 on gut recovery ~9 months)
        fitness.r = rep(1.5, N^2), 
        mu = rep(0.02, N^2),      # "mu", probability of decolonisation (Haggai Bar-Yoseph, JAC, 2016, decreasing colonization rates from 76.7% (95% CI=69.3%–82.8%) at 1 month to 35.2% (95% CI=28.2%–42.9%) at 12 months of follow-up)
        abx.s = rep(0.5, N^2),     # abx.s = amount of s killed by broad spectrum abx s
        abx.r = rep(0, N^2),         # abx.r = amount of r killed by broad spectrum abx r
        p.infect = rep(0.5, N^2),     # "p.infect", probability of being prescribed antibiotics
        p.infect.after = rep(0.02, N^2),    # admission day when cummulative prabability of HAI requiring abx.r is 1
        p.r.day1 = rep(0.5, N^2),    # probability of being prescribed broad spectrum antibiotic on admission 
        p.r.after = rep(0.5, N^2),    # probability of being prescribed broad spectrum antibiotic after admission 
        meanDur = samples[['meanDur']], 
        iterations = rep(iterations, N^2)                         
    )
    
    # Check order of variables - MAKE SURE the variable listing and ORDER MATCHES the variable listing input into run_model
    source('models/get_output_abs_cocarriage5state.R')
    if(all(names(sampled.para)  == formalArgs(run_abs_cocarriage5state))){
        
        old = Sys.time() # get start time
        
        sampled.para.df.list = split(sampled.para, seq(nrow(sampled.para)))
        sampled.para.list = lapply(sampled.para.df.list, unlist)
        
        out = parallel::parLapply(cl, sampled.para.list, function(x) {
            
            print(x)
            
            trialdata = run_abs_cocarriage5state(n.bed= x['n.bed'], max.los = x['max.los'], 
                                                 prop_R= x['prop_R'], prop_r = x['prop_r'], prop_Sr = x['prop_Sr'], prop_S = x['prop_S'], 
                                                 bif = x['bif'], pi_ssr = x['pi_ssr'], repop.s = x['repop.s'], fitness.r = x['fitness.r'], mu = x['mu'], 
                                                 abx.s = x['abx.s'], abx.r = x['abx.r'],
                                                 p.infect = x['p.infect'], p.infect.after = x['p.infect.after'], 
                                                 p.r.day1 = x['p.r.day1'], p.r.after = x['p.r.after'], 
                                                 meanDur = x['meanDur'],
                                                 iterations = x['iterations'])
            
            return(trialdata)
        })

        print(Sys.time() - old) # print elapsed time difference in nice format
    }
    
    
    d = do.call('rbind', out)
    colnames(d) = res.names
    dat[[p]] = cbind(sampled.para, d, rep(p, nrow(d)))
    colnames(dat[[p]])[grep('nrow', colnames(dat[[p]]))] = 'p' 
}


stopCluster(cl)

dat = do.call('rbind', dat)
save(dat, file = paste0("runs/cocarriage5state_transheat_abs_", Sys.Date(), ".Rdata"))
    



