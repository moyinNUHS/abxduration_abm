# get absolute difference between long and short duration
# in terms of new R acquisitions and R carrier prevalence 

source('models/los_abx_matrix_varydur.R')
source('models/summary_los.R')
source('models/model_populationgrowth.R')

get_output_summarystats_populationgrowth <- function(n.bed, max.los, n.day, 
                                                     prop_R, r_thres, K, total_prop,
                                                     pi_ssr, s_growth, fitness.r, r_trans,
                                                     abx.s, abx.r, 
                                                     p.infect, cum.r.1, p.r.day1, p.r.after,
                                                     meanDur, dur.type, timestep, iterations, burn_in){
  
  message(paste0('running ', dur.type, ' duration for ', iterations, ' iterations...'))
  
  # empty matrix to store output - each row is a day, each col is a iteration 
  iter_totalR.thres = matrix(NA, nrow = n.day, ncol = iterations)
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
    colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, total_prop=total_prop, prop_R=prop_R, r_thres=r_thres, K=K)
    
    # update values for each day - list 1 contains S, list 2 contains R
    colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                      pi_ssr=pi_ssr, total_prop=total_prop, K=K, 
                                      fitness.r=fitness.r, r_thres=r_thres,r_trans=r_trans, s_growth=s_growth,
                                      abx.s=abx.s, abx.r=abx.r, timestep=timestep)
    
    # summary of the output - timestep by bed in absolute numbers
    df.R = colo.matrix_filled_iter[[2]] # get matrix containing R 
    
    r_thres_matrix = colo.matrix_filled_iter[[4]]
    
    # for number of people who reached R threshold on a day
    ##   sum of number of people per timestep that reach threshold 
    ##   make a matrix of sum of people per day (days by timestep)
    ##   daily means of number of people who reached R threshold 
    total_overRthres_perday = matrix(rowSums(df.R >= r_thres_matrix), ncol=timestep, byrow=T)
    iter_totalR.thres[, iter] = rowMeans(total_overRthres_perday)
    
    ### New R at discharge when not admitted as R 
    patient.matrix.exlude.burnin = patient.matrix[((burn_in + 1) : n.day), ]   # patient matrix after the first `burn_in` days 
    patient.id.exlude.burnin = unique(as.vector(patient.matrix.exlude.burnin)) # patient ids admitted after the first `burn_in` days 
    admit.positions = c(1, cumsum(los.array[2,]) + 1)[1:length(los.array[2,])]
    admit.positions.exlude.burnin = admit.positions[patient.id.exlude.burnin]  # admission days of patients admitted after the first `burn_in` days 
    discharge.positions = cumsum(los.array[2,])
    discharge.positions.exlude.burnin = discharge.positions[patient.id.exlude.burnin] # discharge days of patients admitted after the first `burn_in` days 
    colo_vector_filled_iter = as.vector(df.R)
    colo_vector_filled_iter_TF = colo_vector_filled_iter >= as.vector(r_thres_matrix) # turn values into T (R carrier) or F (did not reach threshold for R carrier)
    
    if (any(colo_vector_filled_iter_TF == T)) {
      new_R = length(which(colo_vector_filled_iter_TF[admit.positions.exlude.burnin] == F & colo_vector_filled_iter_TF[discharge.positions.exlude.burnin] == T))
      new_R_peradmission = new_R / sum(colo_vector_filled_iter_TF[admit.positions.exlude.burnin] == F)
    } else {
      new_R_peradmission = 0
    }
    iter_newR[iter] = new_R_peradmission
  }
  
  # Discard first `burn_in` days as burn-in
  total_R_periter_perday = rowSums(iter_totalR.thres[(burn_in + 1) : nrow(iter_totalR.thres), , drop = FALSE]) / iterations / n.bed
  totalR = mean(total_R_periter_perday) # mean of the days 
  newR = mean(iter_newR) # mean per iteration
  
  return(c(totalR = totalR, newR = newR))
  
}

# function to run the model 
run_absdiff_populationgrowth <- function(n.bed, max.los, p.infect, cum.r.1, p.r.day1, p.r.after,
                                         K, total_prop,  prop_R, pi_ssr, 
                                         r_trans, fitness.r, r_thres, s_growth,
                                         abx.s, abx.r, short_dur,long_dur){
  
  message('population growth model initiating...')
  
  timestep = 1
  iterations = 100 # from AA tests
  n.day = 300
  burn_in = 150   # from equilibrium graphs 
  
  ############ 
  #### SHORT DURATION 
  
  short_output = get_output_summarystats_populationgrowth(n.bed = n.bed, n.day = n.day, max.los = max.los, cum.r.1 = cum.r.1,
                                                          prop_R = prop_R, r_thres = r_thres, K = K, total_prop = total_prop,
                                                          pi_ssr = pi_ssr,  s_growth = s_growth, fitness.r = fitness.r, r_trans = r_trans,
                                                          abx.r = abx.r, abx.s = abx.s, timestep = timestep, iterations = iterations,
                                                          p.infect = p.infect, p.r.day1 = p.r.day1, p.r.after = p.r.after, 
                                                          meanDur = short_dur, dur.type = 'short', burn_in = burn_in)

  ############ 
  #### LONG DURATION 
  
  long_output = get_output_summarystats_populationgrowth(n.bed = n.bed, n.day = n.day, max.los = max.los, cum.r.1 = cum.r.1,
                                                          prop_R = prop_R, r_thres = r_thres, K = K, total_prop = total_prop,
                                                          pi_ssr = pi_ssr,  s_growth = s_growth, fitness.r = fitness.r, r_trans = r_trans,
                                                          abx.r = abx.r, abx.s = abx.s, timestep = timestep, iterations = iterations,
                                                          p.infect = p.infect, p.r.day1 = p.r.day1, p.r.after = p.r.after, 
                                                          meanDur = long_dur, dur.type = 'long', burn_in = burn_in)
  
  out = c(long_output, short_output, long_output - short_output)
  
  return(out)
}

# names of the output columns
res.names = c('totalR_long', 'newR_long', 
              'totalR_short', 'newR_short', 
              'totalR_diff', 'newR_diff')
