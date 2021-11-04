# get % R carriers amongst those treated  

source('models/los_abx_matrix_varydur.R')
source('models/summary_los.R')
source('models/model_cocarriage5state.R')

run_abs_cocarriage5state <- function(n.bed, max.los, 
                                     prop_R, prop_r, prop_Sr, prop_S,
                                     bif, pi_ssr, repop.s, fitness.r, 
                                     mu, abx.s, abx.r, 
                                     p.infect, p.infect.after, p.r.day1, p.r.after,
                                     meanDur, iterations = 100){
  
  timestep = 1
  n.day = 300
  burn_in = 150    # from equilibrium graphs 
  
  # empty matrix to store output - each row is a day, each col is a iteration 
  
  propstatusperday = list()
  for (c in c('S', 'ss', 'Sr', 'sr', 'sR')){
    
    iter_total = matrix(NA, nrow = n.day, ncol = iterations)
    
    for(iter in 1:iterations){ # for each iteration 
      
      # get matrix of length of stay, abx prescribed, patients admitted
      matrixes = los_abx_table_varydur(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                                       p.infect=p.infect, p.r.day1=p.r.day1, p.r.after = p.r.after, 
                                       p.infect.after=p.infect.after, 
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
      ### per day 
      total_perday = matrix(rowSums(colo.matrix_filled_iter == c), ncol = timestep, byrow=T)
      iter_total[, iter]= rowMeans(total_perday)
      
    }
    
    # Discard first `burn_in` days as burn-in
    total_periter_perday = rowSums(iter_total[(burn_in + 1) : nrow(iter_total), , drop = FALSE]) / iterations / n.bed
    propstatusperday[[c]] = mean(total_periter_perday) # mean of the days 
    
  }
  
  return(propstatusperday)
  
}


# names of the output columns
res.names = paste('prop', c('S', 's', 'Sr', 'sr', 'sR'), 'per day')
