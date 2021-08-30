source('models/los_abx_matrix_withinhost.R')
source('models/summary_los.R')
source('models/model_populationgrowth.R')

run_trialdata_withinhost_populationgrowth <- function(n.bed, n.day, 
                                                     prop_R, r_thres, K, total_prop,
                                                     s_growth, fitness.r, r_trans,
                                                     abx.s, abx.r, short_dur, long_dur, 
                                                     timestep, iterations){
  
  message(paste0('running trialdata_withinhost_populationgrowth for ', iterations, ' iterations...'))
  
  # empty matrix to store output - each row is a day, each col is a iteration 
  obs = array(NA, dim = c(n.day, iterations, 3, 2), 
              dimnames = list(NULL, NULL, c('R', 'S', 'r_thres'), c('short', 'long')))
  
  for(iter in 1:iterations){ # for each iteration 
    
    # get matrix of length of stay, abx prescribed, patients admitted
    matrixes = los_abx_table_withinhost(n.bed = n.bed, n.day = n.day, short_dur = short_dur, long_dur = long_dur)
    patient.matrix = matrixes[[1]]
    abx.matrix.short = matrixes[[2]]
    abx.matrix.long = matrixes[[3]]
    los.array = summary_los(patient.matrix = patient.matrix)
    
    # starting state for all the patients admitted 
    colo.matrix = colo.table(patient.matrix = patient.matrix, los.array=los.array, total_prop=total_prop, prop_R=prop_R, r_thres=r_thres, K=K)
    
    # update values for each day - list 1 contains S, list 2 contains R
    colo.matrix_filled_iter_short = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix.short, colo.matrix=colo.matrix, 
                                            pi_ssr=0, total_prop=total_prop, K=K, 
                                            fitness.r=fitness.r, r_thres=r_thres,r_trans=r_trans, s_growth=s_growth,
                                            abx.s=abx.s, abx.r=abx.r, timestep=timestep)
    
    colo.matrix_filled_iter_long = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix.long, colo.matrix=colo.matrix, 
                                           pi_ssr=0, total_prop=total_prop, K=K, 
                                           fitness.r=fitness.r, r_thres=r_thres,r_trans=r_trans, s_growth=s_growth,
                                           abx.s=abx.s, abx.r=abx.r, timestep=timestep)
    
    
    # summary of the output 
    total_capacity_matrix_short = colo.matrix_filled_iter_short[[3]]
    df.S.short = exp(colo.matrix_filled_iter_short[[1]])/exp(total_capacity_matrix_short) # get matrix containing S 
    df.R.short = exp(colo.matrix_filled_iter_short[[2]])/exp(total_capacity_matrix_short) # get matrix containing R 
    
    total_capacity_matrix_long = colo.matrix_filled_iter_long[[3]]
    df.S.long = exp(colo.matrix_filled_iter_long[[1]])/exp(total_capacity_matrix_long) # get matrix containing S 
    df.R.long = exp(colo.matrix_filled_iter_long[[2]])/exp(total_capacity_matrix_long) # get matrix containing R 
    
    r_thres_prop_long = mean(exp(colo.matrix_filled_iter_long[[4]])/exp(total_capacity_matrix_long)[1,])
    r_thres_prop_short = mean(exp(colo.matrix_filled_iter_short[[4]])/exp(total_capacity_matrix_short)[1,])
    
    out = lapply(list(df.R.short, df.R.long, df.S.short, df.S.long), rowMeans)
    
    obs[, iter, 'R', 'short'] = out[[1]]
    obs[, iter, 'R', 'long'] = out[[2]]
    obs[, iter, 'S', 'short'] = out[[3]]
    obs[, iter, 'S', 'long'] = out[[4]]
    obs[, iter, 'r_thres', 'short'] = rep(r_thres_prop_short, n.day)
    obs[, iter, 'r_thres', 'long'] = rep(r_thres_prop_long, n.day)
    
  }
  
  # get average of all iterations 
  obs.avg = lapply(list(obs[, , 'R', 'short'], obs[, , 'R', 'long'], 
                        obs[, , 'S', 'short'], obs[, , 'S', 'long'],
                        obs[, , 'r_thres', 'short'], obs[, , 'r_thres', 'long']), rowMeans)
  
  obs.avg.df = as.data.frame(obs.avg)
  colnames(obs.avg.df) = c('R_short_perday', 'R_long_perday', 'S_short_perday', 'S_long_perday', 'r_thres_short', 'r_thres_long')
    
  return(obs.avg.df)
  
}
