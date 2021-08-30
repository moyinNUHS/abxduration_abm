# get % R carriers amongst those treated  

source('models/los_abx_matrix_varydur.R')
source('models/summary_los.R')
source('models/model_populationgrowth.R')

run_treated_populationgrowth <- function(n.bed, max.los, 
                                         p.infect, cum.r.1, p.r.day1, p.r.after,
                                         K, total_prop, prop_R, pi_ssr, 
                                         r_trans, fitness.r,
                                         r_thres, s_growth, 
                                         abx.s, abx.r, 
                                         meanDur){
  
  timestep = 1
  iterations = 100 # from AA tests - switch iterations to 1 when doing AA test
  n.day = 300
  burn_in = 150    # from equilibrium graphs 
  
  message(paste0('running population growth model for outcome amongst treated for ', iterations, ' iterations...'))
  
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
    colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, total_prop=total_prop, prop_R=prop_R, r_thres=r_thres, K=K)
    
    # update values for each day - list 1 contains S, list 2 contains R
    colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                      pi_ssr=pi_ssr, total_prop=total_prop, K=K, 
                                      fitness.r=fitness.r, r_thres=r_thres,r_trans=r_trans, s_growth=s_growth,
                                      abx.s=abx.s, abx.r=abx.r, timestep=timestep)
    
    # summary of the output - timestep by bed in absolute numbers
    df.R = colo.matrix_filled_iter[[2]] # get matrix containing R 
    r_thres_matrix = colo.matrix_filled_iter[[4]]
    colo_vector_filled_iter_TF = df.R >= r_thres_matrix # turn values into T (R carrier) or F (did not reach threshold for R carrier)
    
    if (any(colo_vector_filled_iter_TF == T)) {
      patient.matrix.exlude.burnin = patient.matrix[((burn_in + 1) : n.day), ]   # patient matrix after the first `burn_in` days 
      abx.matrix.exlude.burnin = abx.matrix[((burn_in + 1) : n.day), ]                # abx matrix after the first `burn_in` days 
      colo.matrix.exclude.burnin = colo_vector_filled_iter_TF[((burn_in + 1) : n.day), ] # colo matrix after the first `burn_in` days 
      
      abx.list = split(abx.matrix.exlude.burnin, patient.matrix.exlude.burnin)
      abx.count = unlist(lapply(abx.list, function(x){
        x[which(x != 0)] = 1
        unname(unlist(lapply(split(x, cumsum(c(0, diff(x) != 0))), cumsum)))
      }))
      
      for (day in 1:20){ # every day of antibiotic treatment
        iter_propR_treated_delta[day, iter] = mean(as.vector(colo.matrix.exclude.burnin)[which(abx.count == day)] == T) - prop_R
      }
      
      # proportion of R in those treated 
      iter_propR_treated_ward[iter] = sum(colo.matrix.exclude.burnin[which(abx.matrix.exlude.burnin > 0)] == T) / 
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

