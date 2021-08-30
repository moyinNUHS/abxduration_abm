# get % R carriers amongst those treated  

source('models/los_abx_matrix_varydur.R')
source('models/summary_los.R')
source('models/model_cocarriage5state.R')

get_output_summarystats_cocarriage <- function(n.bed, max.los, 
                                               prop_R, prop_r, prop_Sr, prop_S,
                                               bif, pi_ssr, repop.s, fitness.r, 
                                               mu, abx.s, abx.r, 
                                               p.infect, cum.r.1, p.r.day1, p.r.after,
                                               meanDur){
  
  message(paste0('running ', dur.type, ' duration for ', iterations, ' iterations...'))
  
  # empty matrix to store output - each row is a day, each col is a iteration 
  iter_totalR = matrix(NA, nrow = n.day, ncol = iterations)
  iter_newR = c()
  
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
    colo.matrix = colo.table(patient.matrix = patient.matrix, los=los.array, 
                             prop_R = prop_R, prop_r=prop_r, prop_Sr=prop_Sr, prop_S=prop_S)
    
    # update values for each day 
    colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                      pi_ssr=pi_ssr, bif=bif, mu=mu, fitness.r = fitness.r,
                                      repop.s=repop.s, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
    
    # summary of the output 
    ### sRs per day 
    total_R_perday = matrix(rowSums(colo.matrix_filled_iter == "sR"), ncol = timestep, byrow=T)
    iter_totalR[, iter]= rowMeans(total_R_perday)
    
    ### New R at discharge when not admitted as R 
    patient.matrix.exlude.burnin = patient.matrix[((burn_in + 1) : n.day), ]   # patient matrix after the first `burn_in` days 
    patient.id.exlude.burnin = unique(as.vector(patient.matrix.exlude.burnin)) # patient ids admitted after the first `burn_in` days 
    admit.positions = c(1, cumsum(los.array[2,]) + 1)[1:length(los.array[2,])]
    admit.positions.exlude.burnin = admit.positions[patient.id.exlude.burnin]  # admission days of patients admitted after the first `burn_in` days 
    discharge.positions = cumsum(los.array[2,])
    discharge.positions.exlude.burnin = discharge.positions[patient.id.exlude.burnin] # discharge days of patients admitted after the first `burn_in` days 
    colo_vector_filled_iter = as.vector(colo.matrix_filled_iter)
    if ('sR' %in%  colo_vector_filled_iter) {
      new_R = length(which(colo_vector_filled_iter[admit.positions.exlude.burnin] != 'sR' & colo_vector_filled_iter[discharge.positions.exlude.burnin] == 'sR'))
      new_R_peradmission = new_R/ sum(colo_vector_filled_iter[admit.positions.exlude.burnin] != 'sR') # as a proportion of non-carriers admitted 
    } else {
      new_R_peradmission = 0
    }
    
    iter_newR[iter] = new_R_peradmission
    
  }
  
  # Discard first `burn_in` days as burn-in
  total_R_periter_perday = rowSums(iter_totalR[(burn_in + 1) : nrow(iter_totalR), , drop = FALSE]) / iterations / n.bed
  totalR = mean(total_R_periter_perday) # mean of the days 
  
  newR = mean(iter_newR) # mean per iteration
  
  return(c(totalR = totalR, newR = newR))
  
}

# function to run the model 
run_absdiff_cocarriage5state <- function(n.bed, max.los, 
                                         prop_R, prop_r, prop_Sr, prop_S,
                                         bif, pi_ssr, repop.s, fitness.r,
                                         mu, abx.s, abx.r, 
                                         p.infect, cum.r.1, p.r.day1, p.r.after,
                                         short_dur, long_dur){
  
  message('cocarriage 5 state model initiating...')
  
  timestep = 1
  iterations = 100   # from AA tests
  n.day = 300
  burn_in = 150    # from equilibrium graphs 
  
  ############ 
  #### SHORT DURATION 
  
  short_output = get_output_summarystats_cocarriage(n.bed = n.bed, n.day = n.day, max.los = max.los, 
                                                    p.infect = p.infect, p.r.day1 = p.r.day1, p.r.after = p.r.after, 
                                                    cum.r.1 = cum.r.1, 
                                                    meanDur = short_dur, dur.type = 'short',
                                                    prop_R = prop_R, prop_r = prop_r, prop_Sr = prop_Sr, prop_S = prop_S, 
                                                    pi_ssr = pi_ssr, bif = bif, mu = mu, fitness.r = fitness.r,
                                                    repop.s = repop.s, abx.r=abx.r, abx.s=abx.s, timestep = timestep, iterations = iterations, 
                                                    burn_in = burn_in)
  
  ############ 
  #### LONG DURATION 
  long_output = get_output_summarystats_cocarriage(n.bed = n.bed, n.day = n.day, max.los = max.los, 
                                                   p.infect = p.infect, p.r.day1 = p.r.day1, p.r.after = p.r.after, 
                                                   cum.r.1 = cum.r.1, 
                                                   meanDur = long_dur, dur.type = 'long',
                                                   prop_R = prop_R, prop_r = prop_r, prop_Sr = prop_Sr, prop_S = prop_S, 
                                                   pi_ssr = pi_ssr, bif = bif, mu = mu, fitness.r = fitness.r,
                                                   repop.s = repop.s, abx.r=abx.r, abx.s=abx.s, timestep = timestep, iterations = iterations,
                                                   burn_in = burn_in)
  
  
  out = c(long_output, short_output, long_output - short_output)
  
  return(out)
}

# names of the output columns
res.names = c('totalR_long', 'newR_long', 
              'totalR_short', 'newR_short', 
              'totalR_diff', 'newR_diff')
