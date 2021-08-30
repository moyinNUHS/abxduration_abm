# get % R carriers amongst those treated  

source('models/los_abx_matrix_varydur.R')
source('models/summary_los.R')
source('models/model_simple3state.R')

# function to run the model 
run_treated_simple3state <- function(n.bed, max.los, 
                                     prop_R, prop_S,
                                     bif, pi_ssr, repop.s, 
                                     mu, abx.s, abx.r, 
                                     p.infect, cum.r.1, p.r.day1, p.r.after,
                                     meanDur){
  
  timestep = 1
  iterations = 100 # from AA tests - switch iterations to 1 when doing AA test
  n.day = 300
  burn_in = 150    # from equilibrium graphs 
  
  message(paste0('running simple model for outcome amongst treated for ', iterations, ' iterations...'))
  
  # empty matrix to store output - each row is a day, each col is a iteration 
  iter_propR_treated_delta = matrix(NA, nrow = 20, ncol = iterations)
  iter_propR_treated_ward = c()
  
  for(iter in 1:iterations){ # for each iteration 
    
    # get matrix of length of stay, abx prescribed, patients admitted
    matrixes = los_abx_table_varydur(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                                     p.infect=p.infect, p.r.day1=p.r.day1, p.r.after = p.r.after, 
                                     cum.r.1=cum.r.1, 
                                     meanDur = meanDur, timestep=timestep)
    
    
    patient.matrix = matrixes[[1]]
    abx.matrix = matrixes[[2]]
    los.array = summary_los(patient.matrix = patient.matrix)
    
    
    # starting state for all the patients admitted 
    colo.matrix = colo.table(patient.matrix = patient.matrix, los = los.array, 
                             prop_R = prop_R, prop_S = prop_S)
    
    # update values for each day 
    colo.matrix_filled_iter = nextDay(patient.matrix = patient.matrix, los.array = los.array, 
                                      abx.matrix = abx.matrix, colo.matrix = colo.matrix, 
                                      bif = bif, pi_ssr = pi_ssr, repop.s = repop.s, mu = mu, 
                                      abx.s = abx.s, abx.r = abx.r, timestep = timestep)
    
    # summary of the outputs
    ### R carriers amongst those treated with antibiotics
    if ('R' %in%  colo.matrix_filled_iter) {
      
      patient.matrix.exlude.burnin = patient.matrix[((burn_in + 1) : n.day), ]        # patient matrix after the first `burn_in` days 
      abx.matrix.exlude.burnin = abx.matrix[((burn_in + 1) : n.day), ]                # abx matrix after the first `burn_in` days 
      colo.matrix.exclude.burnin = colo.matrix_filled_iter[((burn_in + 1) : n.day), ] # colo matrix after the first `burn_in` days 
      
      abx.list = split(abx.matrix.exlude.burnin, patient.matrix.exlude.burnin)
      abx.count = unlist(lapply(abx.list, function(x){
        x[which(x != 0)] = 1
        unname(unlist(lapply(split(x, cumsum(c(0, diff(x) != 0))), cumsum)))
      }))
   
      for (day in 1:20){ # every day of antibiotic treatment
        iter_propR_treated_delta[day, iter] = mean(as.vector(colo.matrix.exclude.burnin)[which(abx.count == day)] == 'R') - prop_R
      }
      
      # proportion of R in those treated 
      iter_propR_treated_ward[iter] = sum(colo.matrix.exclude.burnin[which(abx.matrix.exlude.burnin > 0)] == 'R') / 
        sum(abx.matrix.exlude.burnin > 0)
      
    } 

  }
  
  propR_treated_delta = rowMeans(iter_propR_treated_delta, na.rm = T) # mean per iteration
  propR_treated_delta_diffavg = mean(diff(propR_treated_delta[!is.na(propR_treated_delta)]))
  propR_treated_ward = mean(iter_propR_treated_ward)
  
  return(c(propR_treated_delta, propR_treated_delta_diffavg, propR_treated_ward))
  
}

# names of the output columns
res.names = c(paste0('propR_treated_delta', 1:20), 'propR_treated_delta_diffavg', 'propR_treated_ward')
