#########################################################################
#######Effect of antibiotic duration on hospitalised patients############
######################Show difference over time##########################
#########################################################################

default.values <- function(model, scenario){
  
  #taken from source('plots/default_params.R')
  common = t(as.matrix(common.para))
  
  specific = t(as.matrix(specific.para[[model]]))
  
  valuestorepeat = cbind(common, specific)
  if (scenario == 'cre') {valuestorepeat[,'abx.r'] = 0}
  
  #empty matrix 
  m = matrix(rep(unlist(valuestorepeat), iterations), byrow = T,
             ncol = length(valuestorepeat), nrow = iterations, 
             dimnames = list(NULL, colnames(valuestorepeat)))
  
  return(m)
  
}


get.stability.data <- function (model, iterations, scenario) { 
  
  ##get default values to run models 
  feed = default.values(model, scenario)
  
  ##source model code 
  source('models/los_abx_matrix_varydur.R')
  source('models/summary_los.R')
  source(paste0('models/model_', model, '.R'))
  
  ## run models 
  if(model == "simple3state"){
    
    out = apply(feed, 1, function(x){
      
      ############ 
      #### SHORT DURATION 
      
      # get matrix of length of stay, abx prescribed, patients admitted
      matrixes = los_abx_table_varydur(n.bed = x[['n.bed']], n.day = x[['n.day']], max.los = x[['max.los']], 
                                       p.infect = x[['p.infect']], p.r.day1 = x[['p.r.day1']], p.r.after = x[['p.r.after']], 
                                       cum.r.1 = x[['cum.r.1']], 
                                       meanDur = x[['short_dur']], timestep = x[['timestep']])
      patient.matrix = matrixes[[1]]
      abx.matrix = matrixes[[2]]
      los.array = summary_los(patient.matrix = patient.matrix)
      
      # starting state for all the patients admitted 
      colo.matrix = colo.table(patient.matrix = patient.matrix, los = los.array, 
                               prop_R = x[['prop_R']], prop_S = x[['prop_S']])
      
      # update values for each day 
      colo_table_filled_iter = nextDay(patient.matrix = patient.matrix, los.array = los.array, 
                                       abx.matrix = abx.matrix, colo.matrix = colo.matrix, 
                                       bif = x[['bif']], pi_ssr = x[['pi_ssr']], repop.s = x[['repop.s']], mu = x[['mu']], 
                                       abx.s = x[['abx.s']], abx.r = x[['abx.r']], timestep = x[['timestep']])
      
      # summary of the output 
      ### Rs per day 
      total_R_perday = matrix(rowSums(colo_table_filled_iter == "R"), ncol = x[['timestep']], byrow=T)
      iter_totalR = rowMeans(total_R_perday) # average for each time step
      
      cumsum_totalR_short = cumsum(iter_totalR) # sum over columns (days)
      
      ############ 
      #### LONG DURATION 
      
      # get matrix of length of stay, abx prescribed, patients admitted
      matrixes = los_abx_table_varydur(n.bed = x[['n.bed']], n.day = x[['n.day']], max.los = x[['max.los']], 
                                       p.infect = x[['p.infect']], p.r.day = x[['p.r.day1']], p.r.after = x[['p.r.after']], 
                                       cum.r.1 = x[['cum.r.1']], 
                                       meanDur = x[['long_dur']], timestep = x[['timestep']])
      patient.matrix = matrixes[[1]]
      abx.matrix = matrixes[[2]]
      los.array = summary_los(patient.matrix = patient.matrix)
      
      # starting state for all the patients admitted 
      colo.matrix = colo.table(patient.matrix = patient.matrix, los = los.array, 
                               prop_R = x[['prop_R']], prop_S = x[['prop_S']])
      
      # update values for each day 
      colo_table_filled_iter = nextDay(patient.matrix = patient.matrix, los.array = los.array, 
                                       abx.matrix = abx.matrix, colo.matrix = colo.matrix, 
                                       bif = x[['bif']], pi_ssr = x[['pi_ssr']], repop.s = x[['repop.s']], mu = x[['mu']], 
                                       abx.s = x[['abx.s']], abx.r = x[['abx.r']], timestep = x[['timestep']])
      
      # summary of the output 
      ### Rs per day 
      total_R_perday = matrix(rowSums(colo_table_filled_iter == "R"), ncol = x[['timestep']], byrow=T)
      iter_totalR = rowMeans(total_R_perday) # average for each time step
      
      cumsum_totalR_long = cumsum(iter_totalR)
      
      return(cumsum_totalR_long - cumsum_totalR_short)
      
    })
    
  } else if (model == "cocarriage5state") {
    
    out =  apply(feed, 1, function(x){
      
      ############ 
      #### SHORT DURATION 
      
      # get matrix of length of stay, abx prescribed, patients admitted
      matrixes = los_abx_table_varydur(n.bed= x[['n.bed']], n.day= x[['n.day']], max.los= x[['max.los']], 
                                       p.infect= x[['p.infect']], p.r.day1= x[['p.r.day1']], p.r.after = x[['p.r.after']], 
                                       cum.r.1= x[['cum.r.1']], 
                                       meanDur= x[['short_dur']], timestep= x[['timestep']])
      patient.matrix=matrixes[[1]]
      abx.matrix = matrixes[[2]]
      los.array = summary_los(patient.matrix = patient.matrix)
      
      # starting state for all the patients admitted 
      colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                               prop_R= x[['prop_R']], prop_r= x[['prop_r']], prop_Sr= x[['prop_Sr']], prop_S= x[['prop_S']])
      
      # update values for each day 
      colo.matrix_filled_iter = nextDay(patient.matrix = patient.matrix, abx.matrix = abx.matrix, colo.matrix = colo.matrix, 
                                        pi_ssr= x[['pi_ssr']], bif= x[['bif']], mu= x[['mu']], fitness.r = x[['fitness.r']],
                                        repop.s= x[['repop.s']], abx.r= x[['abx.r']], abx.s= x[['abx.s']], timestep= x[['timestep']])
      
      # summary of the output 
      ### sRs per day 
      total_sR_perday = matrix(rowSums(colo.matrix_filled_iter == "sR"), ncol = x[['timestep']], byrow=T)
      iter_totalsR = rowMeans(total_sR_perday)
      
      cumsum_totalR_short = cumsum(iter_totalsR)
      
      ############ 
      #### LONG DURATION 
      
      # get matrix of length of stay, abx prescribed, patients admitted
      matrixes = los_abx_table_varydur(n.bed= x[['n.bed']], n.day= x[['n.day']], max.los= x[['max.los']], 
                                       p.infect= x[['p.infect']], p.r.day1= x[['p.r.day1']], p.r.after = x[['p.r.after']], 
                                       cum.r.1= x[['cum.r.1']], 
                                       meanDur= x[['long_dur']], timestep= x[['timestep']])
      patient.matrix = matrixes[[1]]
      abx.matrix = matrixes[[2]]
      los.array = summary_los(patient.matrix=patient.matrix)
      
      # starting state for all the patients admitted 
      colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                               prop_R= x[['prop_R']], prop_r= x[['prop_r']], prop_Sr= x[['prop_Sr']], prop_S= x[['prop_S']])
      
      # update values for each day 
      colo.matrix_filled_iter = nextDay(patient.matrix= patient.matrix, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                        pi_ssr= x[['pi_ssr']], bif= x[['bif']], mu= x[['mu']], fitness.r = x[['fitness.r']],
                                        repop.s= x[['repop.s']], abx.r= x[['abx.r']], abx.s= x[['abx.s']], timestep= x[['timestep']])
      
      #summary of the output 
      ### sRs per day 
      total_sR_perday = matrix(rowSums(colo.matrix_filled_iter == "sR"), ncol = x[['timestep']], byrow=T)
      iter_totalsR = rowMeans(total_sR_perday) # average for each time step
      
      cumsum_totalR_long = cumsum(iter_totalsR)
      
      return(cumsum_totalR_long - cumsum_totalR_short)
      
    })
    
  } else { #population growth model 
    
    out =  apply(feed, 1, function(x){
      
      # get matrix of length of stay, abx prescribed, patients admitted
      matrixes = los_abx_table_varydur(n.bed= x[['n.bed']], n.day= x[['n.day']], max.los= x[['max.los']], 
                                       p.infect= x[['p.infect']], p.r.day1= x[['p.r.day1']], p.r.after = x[['p.r.after']],
                                       cum.r.1= x[['cum.r.1']], 
                                       meanDur= x[['short_dur']], timestep= x[['timestep']])
      patient.matrix=matrixes[[1]]
      abx.matrix=matrixes[[2]]
      los.array = summary_los(patient.matrix=patient.matrix)
      
      # starting state for all the patients admitted 
      colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, total_prop= x[['total_prop']], prop_R= x[['prop_R']], r_thres= x[['r_thres']], K= x[['K']])
      
      # update values for each day - list 1 contains S, list 2 contains R
      colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                        pi_ssr= x[['pi_ssr']], total_prop = x[['total_prop']], K= x[['K']], 
                                        fitness.r= x[['fitness.r']], r_thres = x[['r_thres']],r_trans= x[['r_trans']], s_growth= x[['s_growth']],
                                        abx.s= x[['abx.s']], abx.r= x[['abx.r']], timestep= x[['timestep']])
      
      # summary of the output - timestep by bed in absolute numbers
      df.R = colo.matrix_filled_iter[[2]] # get matrix containing R 
      
      r_thres_mat = colo.matrix_filled_iter[[4]]
      
      # for number of people who reached R threshold on a day
      ##   sum of number of people per timestep that reach threshold 
      ##   make a matrix of sum of people per day (days by timestep)
      ##   daily means of number of people who reached R threshold 
      total_overRthres_perday = matrix(rowSums(df.R >= r_thres_mat), ncol = x[['timestep']], byrow=T)
      iter_totalR.thres = rowMeans(total_overRthres_perday)
      
      
      cumsum_totalR_short = cumsum(iter_totalR.thres)
      
      ############ 
      #### LONG DURATION 
      
      # get matrix of length of stay, abx prescribed, patients admitted
      matrixes = los_abx_table_varydur(n.bed= x[['n.bed']], n.day= x[['n.day']], max.los= x[['max.los']], 
                                       p.infect= x[['p.infect']], p.r.day1= x[['p.r.day1']], p.r.after = x[['p.r.after']],
                                       cum.r.1= x[['cum.r.1']], 
                                       meanDur= x[['long_dur']], timestep= x[['timestep']])
      patient.matrix=matrixes[[1]]
      abx.matrix=matrixes[[2]]
      los.array = summary_los(patient.matrix=patient.matrix)
      
      # starting state for all the patients admitted
      colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, total_prop= x[['total_prop']], prop_R= x[['prop_R']], r_thres= x[['r_thres']], K= x[['K']])
      
      # update values for each day 
      colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                        pi_ssr= x[['pi_ssr']], total_prop= x[['total_prop']], K= x[['K']], r_trans= x[['r_trans']],
                                        fitness.r= x[['fitness.r']], r_thres= x[['r_thres']], s_growth= x[['s_growth']],
                                        abx.s= x[['abx.s']], abx.r= x[['abx.r']], timestep= x[['timestep']])
      
      # summary of the output - timestep by bed in absolute numbers
      df.R = colo.matrix_filled_iter[[2]] # get matrix containing R 
      
      r_thres_mat = colo.matrix_filled_iter[[4]]
      
      # for number of people who reached R threshold on a day
      ##   sum of number of people per timestep that reach threshold 
      ##   make a matrix of sum of people per day (days by timestep)
      ##   daily means of number of people who reached R threshold 
      total_overRthres_perday = matrix(rowSums(df.R >= r_thres_mat), ncol=x[['timestep']], byrow=T)
      iter_totalR.thres = rowMeans(total_overRthres_perday)
      
      cumsum_totalR_long = cumsum(iter_totalR.thres)
      
      return(cumsum_totalR_long - cumsum_totalR_short)
      
    })
    
  }
  
  per.day = apply(out, 2, function(x) {x / 1:nrow(out)})
  
  stabilitydata = data.frame(x = 1:(nrow(out)), y = per.day) # divide the cum
  colnames(stabilitydata) = c('x', 1:(ncol(stabilitydata) - 1))
  stabilitydata.melt = melt(stabilitydata, id.vars = 'x', variable_name = 'iter')
  stabilitydata.melt$model = model 
  stabilitydata.melt$scenario = scenario
  
  return(stabilitydata.melt)
  
}



