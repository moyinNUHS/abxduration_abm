################################################################################################
######### Effect of antibiotic duration on resistance carriage in hospitalised patients#########
####### functions to run models for magnitude of change in resistance prevalence ###############
################################################################################################

require(pse) #load pse package for Latin Hypercube
require(sensitivity) #load sensitivity package for sensitivity analysis 
require(parallel) # load parallel processing package to use multiple cores on computer (or cluster)

#fixed parameters 
n.bed.min = 15 
n.bed.max = 15.0001
max.los.min = 3 
max.los.max = 7 
short_dur.min = 0
short_dur.max = 0.00001
long_dur.min = 1
long_dur.max = 21

run_magnitude_simple <- function(pi_ssr.min, pi_ssr.max, 
                                 abx.r.min, abx.r.max, 
                                 abx.s.min, abx.s.max, 
                                 p.infect.min, p.infect.max,          
                                 cum.r.1.min, cum.r.1.max,
                                 prop_R.min, prop_R.max,
                                 p.r.day1.min, p.r.day1.max) {
  
  n.bed.min = 15; n.bed.max = 15.0001
  max.los.min = 3; max.los.max = 7 
  short_dur.min = 0; short_dur.max = 0.0001
  long_dur.min = 1; long_dur.max = 21
  
  cl <- makeCluster(12, outfile = paste0('../../../../Desktop/parallel_error_simple_', Sys.time(), '.txt'))
  
  model <- 'simple'
  
  clusterCall(cl, function() {source('model_simple.R')}) # source functions on all cores
  
  modelRun.simple <- function (data.df) { #data.df is a dataframe of the parameter values in columns 
    return(mapply(diff_prevalence, 
                  data.df[,1], data.df[,2], data.df[,3], 
                  data.df[,4], data.df[,5], data.df[,6], 
                  data.df[,7], data.df[,8], data.df[,9], 
                  data.df[,10], data.df[,11], data.df[,12], 
                  data.df[,13],  data.df[,14], data.df[,15]
    ))
  }
  
  ################################## Define parameters and run LHS ####################################
  
  parameters <- list(
    c("qunif", list(min=n.bed.min, max=n.bed.max), "n.bed"),              #"n.bed", number of beds in the ward
    c("qunif", list(min=max.los.min, max=max.los.max), "max.los"),       #"mean.max.los", mean of length of stay
    c("qunif", list(min=prop_R.min, max=prop_R.max), "prop_R"),    #"prob_StartBact_R",probability of initial carriage of resistant organisms
    c("qunif", list(min=0, max=1), "prop_S"),         #"prop_S", proportion of S in the population of S and ss
    c("qunif", list(min=0, max=1), "bif"),                 #"bif", bacterial interference factor
    c("qunif", list(min=pi_ssr.min, max=pi_ssr.max), "pi_ssr"),            # "pi_ssr" probability of being transmitted r to ss (ss—> ssr)
    c("qunif", list(min=0.02, max=0.12), "repop.s"),     # "repop.s1" probability of ss repopulated to S (Palleja, Nature Biology, 2018 on gut recovery ~9 months)
    c("qunif", list(min=0.002, max=0.02), "mu"),         # "mu_r", probability of decolonisation (Haggai Bar-Yoseph, JAC, 2016, decreasing colonization rates from 76.7% (95% CI=69.3%–82.8%) at 1 month to 35.2% (95% CI=28.2%–42.9%) at 12 months of follow-up)
    c("qunif", list(min=abx.s.min, max=abx.s.max), "abx.s"),           # "abx.s", probability of S becoming ss after being on narrow spectrum antibiotics
    c("qunif", list(min=abx.r.min, max=abx.r.max), "abx.r"),           # "abx.r", probability of R becoming ss after being on broad spectrum antibiotics
    c("qunif", list(min=p.infect.min, max=p.infect.max), "p.infect"),          # "p.infect", probability of being prescribed antibiotics
    c("qunif", list(min=cum.r.1.min, max=cum.r.1.max), "cum.r.1"),        # admission day when cummulative prabability of HAI requiring abx.r is 1
    c("qunif", list(min=p.r.day1.min, max= p.r.day1.max), "p.r.day1"),          # probability of being prescribed broad spectrum antibiotic on admission 
    c("qunif", list(min=short_dur.min, max=short_dur.max), "short_dur"),           # "short_dur", mean short duration of antibiotics (normal distribution)
    c("qunif", list(min=long_dur.min, max=long_dur.max), "long_dur")           # "long_dur", mean long duration of antibiotics (normal distribution)
  )  
  
  q <- unlist(lapply(parameters, function(l) l[[1]]))
  q.arg <- lapply(parameters, function(l) l[2:3])
  factors <- unlist(lapply(parameters, function(l) l[[4]]))
  
  # Use the LHD function to generate a hypercube 
  N = 200
  LHS.simple = LHS(modelRun.simple, factors, N=N, q, q.arg, nboot=100, cl=cl) #N is the size of the hypercube
  image_name = paste0("magnitude_", model, "_", N, "_abxr", abx.r.min, "_pissr", pi_ssr.min, "_pinfect", p.infect.min, format(Sys.time(), "%d%b%Y_%H%M%Z"))
  save(LHS.simple, file = paste0("./runs/", image_name, ".Rdata"))
  
}

run_magnitude_binary <- function(pi_ssr.min, pi_ssr.max, 
                                 abx.r.min, abx.r.max, 
                                 abx.s.min, abx.s.max, 
                                 p.infect.min, p.infect.max,          
                                 cum.r.1.min, cum.r.1.max,
                                 prop_R.min, prop_R.max,
                                 p.r.day1.min, p.r.day1.max) {
  
  n.bed.min = 15; n.bed.max = 15.0001
  max.los.min = 3; max.los.max = 7 
  short_dur.min = 0; short_dur.max = 0.0001
  long_dur.min = 1; long_dur.max = 21
  
  cl <-makeCluster(12, outfile=paste0('/Users/moyin/Desktop/parallel_error_binary_', Sys.time(), '.txt'))
  
  model <- 'binary'
  
  clusterCall(cl, function() {source('model_binary.R')})
  
  modelRun.binary <- function (data.df) { #data.df is a dataframe of the parameter values in columns 
    return(mapply(diff_prevalence, 
                  data.df[,1], data.df[,2], data.df[,3], data.df[,4], data.df[,5], 
                  data.df[,6], data.df[,7], data.df[,8], data.df[,9], 
                  data.df[,10], data.df[,11], data.df[,12], data.df[,13], data.df[,14], data.df[,15], 
                  data.df[,16], data.df[,17], data.df[,18]
    ))
  }
  
  ################################## Define parameters and run LHS ####################################
  parameters <- list(
    c("qunif", list(min=n.bed.min, max=n.bed.max), "n.bed"),              #"n.bed", number of beds in the ward
    c("qunif", list(min=max.los.min, max=max.los.max), "max.los"),       #"mean.max.los", mean of length of stay
    c("qunif", list(min=prop_R.min, max=prop_R.max), "prop_R"),    #"prob_StartBact_R",probability of initial carriage of resistant organisms
    c("qunif", list(min=0, max=1), "prop_r"),         #proportion of S in (S+s): prob_start_S <- prop_S_nonR*(1-prob_R)
    c("qunif", list(min=0, max=1), "prop_Sr"),        #proportion of Sr in (r+R): prob_start_Sr <- prop_Sr_inR*prob_R
    c("qunif", list(min=0, max=1), "prop_S"),         #proportion of sr in (r+r): prob_start_sr <- prop_sr_inR*prob_R
    c("qunif", list(min=0, max=1), "bif"),            #bacterial interference factor (pi_ssr = pi_r1 * bif )
    c("qunif", list(min=pi_ssr.min, max=pi_ssr.max), "pi_ssr"),            # "pi_ssr" probability of being transmitted r to ss (ss—> ssr)
    c("qunif", list(min=0.02, max=0.12), "repop.s"),  #probability of regrowth of S  (s—>S)
    c("qunif", list(min=0.02, max=0.15), "repop.r"),  #probability of regrowth of s (sr—> sR)
    c("qunif", list(min=0.002, max=0.02), "mu"),      #probability of being decolonised to S (Sr—> S) 
    c("qunif", list(min=abx.s.min, max=abx.s.max), "abx.s"),           # "abx.s", probability of S becoming ss after being on narrow spectrum antibiotics
    c("qunif", list(min=abx.r.min, max=abx.r.max), "abx.r"),           # "abx.r", probability of R becoming ss after being on broad spectrum antibiotics
    c("qunif", list(min=p.infect.min, max=p.infect.max), "p.infect"),          # "p.infect", probability of being prescribed antibiotics
    c("qunif", list(min=cum.r.1.min, max=cum.r.1.max), "cum.r.1"),        # admission day when cummulative prabability of HAI requiring abx.r is 1
    c("qunif", list(min=p.r.day1.min, max= p.r.day1.max), "p.r.day1"),          # probability of being prescribed broad spectrum antibiotic on admission 
    c("qunif", list(min=short_dur.min, max=short_dur.max), "short_dur"),           # "short_dur", mean short duration of antibiotics (normal distribution)
    c("qunif", list(min=long_dur.min, max=long_dur.max), "long_dur")           # "long_dur", mean long duration of antibiotics (normal distribution)
  )
  
  q <- unlist(lapply(parameters, function(l) l[[1]]))
  q.arg <- lapply(parameters, function(l) l[2:3])
  factors <- unlist(lapply(parameters, function(l) l[[4]]))
  
  # Use the LHD function to generate a hypercube 
  N = 200
  LHS.binary <- LHS(modelRun.binary, factors, N=N, q, q.arg, nboot=100, cl=cl)
  image_name = paste0("magnitude_", model, "_", N, "_abxr", abx.r.min, "_pissr", pi_ssr.min, "_pinfect", p.infect.min, format(Sys.time(), "%d%b%Y_%H%M%Z"))
  save(LHS.binary, file = paste0("./runs/", image_name, ".Rdata"))
  
}

run_magnitude_frequency <- function(pi_ssr.min, pi_ssr.max, 
                                    abx.r.min, abx.r.max, 
                                    abx.s.min, abx.s.max, 
                                    p.infect.min, p.infect.max,          
                                    cum.r.1.min, cum.r.1.max,
                                    prop_R.min, prop_R.max,
                                    p.r.day1.min, p.r.day1.max) {
  
  n.bed.min = 15; n.bed.max = 15.0001
  max.los.min = 3; max.los.max = 7 
  short_dur.min = 0; short_dur.max = 0.0001
  long_dur.min = 1; long_dur.max = 21
  
  cl <-makeCluster(12, outfile=paste0('/Users/moyin/Desktop/parallel_error_frequency_', Sys.time(), '.txt'))
  
  model <- 'frequency'
  # source functions on all cores
  clusterCall(cl, function() {source('model_frequency.R')})
  
  modelRun.freq <- function (data.df) { #data.df is a dataframe of the parameter values in columns 
    return(mapply(diff_prevalence, 
                  data.df[,1], data.df[,2], data.df[,3], 
                  data.df[,4], data.df[,5], data.df[,6], 
                  data.df[,7], data.df[,8], data.df[,9], 
                  data.df[,10], data.df[,11], data.df[,12], 
                  data.df[,13], data.df[,14], data.df[,15], 
                  data.df[,16], data.df[,17]
    ))
  }
  
  ################################## Define parameters and run LHS ####################################
  parameters <- list(c("qunif", list(min=n.bed.min, max=n.bed.max), "n.bed"),              #"n.bed", number of beds in the ward
                     c("qunif", list(min=max.los.min, max=max.los.max), "max.los"),        #"mean.max.los", mean of length of stay
                     c("qunif", list(min=p.infect.min, max=p.infect.max), "p.infect"),     # "p.infect", probability of being prescribed antibiotics
                     c("qunif", list(min=cum.r.1.min, max=cum.r.1.max), "cum.r.1"),        # admission day when cummulative prabability of HAI requiring abx.r is 1
                     c("qunif", list(min=p.r.day1.min, max= p.r.day1.max), "p.r.day1"),    # probability of being prescribed broad spectrum antibiotic on admission 
                     c("qunif", list(min=18, max=24), "K"),                                # gut holding capacity, on log scale, largest R number possible is exp(300) - typical colonic bacteria 10^14 number/mL content https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4991899/
                     c("qunif", list(min=0.1, max=0.9), "total_prop"),                     # total proportion of holding capacity, the starting amount of enterobacteriaceae
                     c("qunif", list(min=prop_R.min, max=prop_R.max), "prop_R"),           #"prob_R",probability of initial carriage of resistant organisms
                     c("qunif", list(min=pi_ssr.min, max=pi_ssr.max), "pi_ssr"),           # "pi_ssr" probability of being transmitted r to ss (ss—> ssr)
                     c("qunif", list(min=2, max=7), "r_trans"),                            # r_trans = mean amunt of R transmitted
                     c("qunif", list(min=0.3, max=1.4), "r_growth"),                       # r_growth = growth constant for logistic growth
                     c("qunif", list(min=8, max=13), "r_thres"),                           # r_thres = threshold amount of bacteria before R can be transmitted
                     c("qunif", list(min=0.03,max=0.3), "s_growth"),                       # s_growth = amount transmitted on log scale
                     c("qunif", list(min=abx.s.min, max=abx.s.max), "abx.s"),              # "abx.s", probability of S becoming ss after being on narrow spectrum antibiotics
                     c("qunif", list(min=abx.r.min, max=abx.r.max), "abx.r"),              # "abx.r", probability of R becoming ss after being on broad spectrum antibiotics
                     c("qunif", list(min=short_dur.min, max=short_dur.max), "short_dur"),  # "short_dur", mean short duration of antibiotics (normal distribution)
                     c("qunif", list(min=long_dur.min, max=long_dur.max), "long_dur")      # "long_dur", mean long duration of antibiotics (normal distribution)
  )
  
  q <- unlist(lapply(parameters, function(l) l[[1]]))
  q.arg <- lapply(parameters, function(l) l[2:3])
  factors <- unlist(lapply(parameters, function(l) l[[4]]))
  
  # Use the LHD function to generate a hypercube 
  N = 200
  LHS.freq = LHS(modelRun.freq, factors, N=N, q, q.arg, nboot=100, cl=cl)
  image_name = paste0("magnitude_", model, "_", N, "_abxr", abx.r.min, "_pissr", pi_ssr.min, "_pinfect", p.infect.min, format(Sys.time(), "%d%b%Y_%H%M%Z"))
  save(LHS.freq, file = paste0("./runs/", image_name, ".Rdata"))
  
}



