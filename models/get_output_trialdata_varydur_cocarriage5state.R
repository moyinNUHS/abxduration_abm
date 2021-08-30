source('models/los_abx_matrix_varydur.R')
source('models/summary_los.R')
source('models/model_cocarriage5state.R')

# get data per iteration in terms of patient matrix, los, abx and daily states 
get_output_trialdata_varydur_cocarriage <- function(n.bed, max.los, n.day, 
                                                    prop_R, prop_r, prop_Sr, prop_S,
                                                    bif, pi_ssr, repop.s, fitness.r,
                                                    mu, abx.s, abx.r, 
                                                    p.infect, cum.r.1, p.r.day1, p.r.after,
                                                    meanDur, dur.type, timestep, iterations, burn_in){
  
  message(paste0('running ', dur.type, ' duration for ', iterations, ' iterations...'))
  
  # set up empty matrix to record outcomes 
  outcomes.to.monitor.perday = c('prop_R_perday', 'prop_newR_perday') # per day during los
  outcomes.to.monitor.perpatient = c('abxdur.perpatient', 'abx1dur.perpatient', 'abx2dur.perpatient')
  out.perday = array(NA, dim = c(30, iterations, length(outcomes.to.monitor.perday)), # present first 30 days of stay for each patient
                     dimnames = list(NULL, NULL, outcomes.to.monitor.perday))
  out.perpatient = array(NA, dim = c(iterations, length(outcomes.to.monitor.perpatient)), 
                         dimnames = list(NULL, outcomes.to.monitor.perpatient))
  
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
                                      pi_ssr=pi_ssr, bif=bif, mu=mu, fitness.r=fitness.r, 
                                      repop.s=repop.s, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
    
    
    # remove burn-in period 
    days.after.burn.in = (burn_in + 1) : n.day
    patient.matrix.after.burn.in = patient.matrix[days.after.burn.in,]
    pt.id.after.burn.in = unique(as.vector(patient.matrix.after.burn.in))
    abx.matrix.after.burn.in = abx.matrix[days.after.burn.in,]
    carriage.status.after.burn.in = colo.matrix_filled_iter[days.after.burn.in,]
    
    # get each patient's stay in lists 
    los.after.burn.in = rle(as.vector(patient.matrix.after.burn.in))$lengths # duration of los 
    abx.per.patient.after.burn.in = split(as.vector(abx.matrix.after.burn.in), rep.int(seq_along(los.after.burn.in), los.after.burn.in))
    carriage.status.per.patient.after.burn.in = split(as.vector(carriage.status.after.burn.in), rep.int(seq_along(los.after.burn.in), los.after.burn.in))
    carriage.status.per.patient.after.burn.in.samelength = lapply(carriage.status.per.patient.after.burn.in, `length<-`, 30)
    carriage.status.per.patient.df = do.call('rbind', carriage.status.per.patient.after.burn.in.samelength) # fill NA till same length of stay 
    carriage.status.per.patient.df.cut = carriage.status.per.patient.df[, 1:30] 
    
    ####################
    # get summary data #
    ####################
    
    # PATIENTS 
    ### proportion of R carriers for each day during los of 30 days
    out.perday[ , iter, 'prop_R_perday'] = colMeans(carriage.status.per.patient.df.cut == 'sR', na.rm = T)
    
    ### proportion of new acquisitions each day during los of 30 days
    new_R_day_per_patient = lapply(carriage.status.per.patient.after.burn.in, function(status.per.day){
      R_days = c(0, which(status.per.day == 'sR'))
      acquisition_days = R_days[(which(diff(R_days) > 1) + 1)]
      acquisition_days_less30 = acquisition_days[which(acquisition_days < 30)]
      
      stay = rep(0, 30)
      stay[acquisition_days_less30] = 1
      
      return(stay)
    })
    new_R_matrix = do.call('rbind', new_R_day_per_patient)
    R_acquisitons_ontheday = colSums(new_R_matrix == 1, na.rm = T)[-1]
    nonR_onthepreviousday = colSums(new_R_matrix == 0, na.rm = T) [-30]
    out.perday[ , iter, 'prop_newR_perday'] = c(0, R_acquisitons_ontheday / nonR_onthepreviousday)
    
    # ABX USED 
    ### abx duration for each patient 
    abx.dur = lapply(abx.per.patient.after.burn.in, function(abx.per.patient){
      c(sum(abx.per.patient == 1), sum(abx.per.patient == 2), sum(abx.per.patient > 0))
    })
    abx.dur.df = do.call('rbind', abx.dur)
    out.perpatient[iter, 'abx1dur.perpatient'] = mean(abx.dur.df[,1], na.rm = T)
    out.perpatient[iter, 'abx2dur.perpatient'] = mean(abx.dur.df[,2], na.rm = T)
    out.perpatient[iter, 'abxdur.perpatient'] = mean(abx.dur.df[,3], na.rm = T)
    
  }
  
  out = list(prop_R_perday_perparameterset = rowMeans(out.perday[, , 'prop_R_perday'], na.rm = T), 
             prop_newR_perday_perparameterset = rowMeans(out.perday[, , 'prop_newR_perday'], na.rm = T),
             abxdur_perpatient_perparameterset = mean(out.perpatient[,'abxdur.perpatient'], na.rm = T),
             abx1dur_perpatient_perparameterset = mean(out.perpatient[,'abx1dur.perpatient'], na.rm = T),
             abx2dur_perpatient_perparameterset = mean(out.perpatient[,'abx2dur.perpatient'], na.rm = T))
  
  
  return(out)
  
}


#################### run the model using the above functions  #####################
run_trialdata_varydur_cocarriage5state <- function(n.bed, max.los, 
                                                   prop_R, prop_r, prop_Sr, prop_S,
                                                   bif, pi_ssr, repop.s, fitness.r,
                                                   mu, abx.s, abx.r, 
                                                   p.infect, cum.r.1, p.r.day1, p.r.after,
                                                   short_dur, long_dur){
  
  message('cocarriage 5 state model initiating...')
  
  timestep = 1
  iterations = 100 # from AA tests
  n.day = 300
  burn_in = 150    # from equilibrium graphs 
  
  ############ 
  #### SHORT DURATION 
  
  short_output = get_output_trialdata_varydur_cocarriage(n.bed = n.bed, n.day = n.day, max.los = max.los, 
                                                         p.infect = p.infect, p.r.day1 = p.r.day1, p.r.after = p.r.after, 
                                                         cum.r.1 = cum.r.1, 
                                                         meanDur = short_dur, dur.type = 'short',
                                                         prop_R = prop_R, prop_r = prop_r, prop_Sr = prop_Sr, prop_S = prop_S, 
                                                         pi_ssr=pi_ssr, bif=bif, mu=mu, fitness.r=fitness.r, 
                                                         repop.s=repop.s, abx.r=abx.r, abx.s=abx.s, timestep=timestep, iterations = iterations, 
                                                         burn_in = burn_in)
  
  ############ 
  #### LONG DURATION 
  long_output = get_output_trialdata_varydur_cocarriage(n.bed = n.bed, n.day = n.day, max.los = max.los, 
                                                        p.infect = p.infect, p.r.day1 = p.r.day1, p.r.after = p.r.after, 
                                                        cum.r.1 = cum.r.1, 
                                                        meanDur = long_dur, dur.type = 'long',
                                                        prop_R = prop_R, prop_r = prop_r, prop_Sr = prop_Sr, prop_S = prop_S, 
                                                        pi_ssr = pi_ssr, bif = bif, mu = mu, fitness.r = fitness.r, 
                                                        repop.s = repop.s, abx.r = abx.r, abx.s = abx.s, timestep = timestep, iterations = iterations, 
                                                        burn_in = burn_in)
  
  out = list(para = c(n.bed = n.bed, max.los = max.los, 
                      prop_R = prop_R, prop_r = prop_r, prop_Sr = prop_Sr, prop_S = prop_S,
                      bif = bif, pi_ssr = pi_ssr, repop.s = repop.s, fitness.r = fitness.r,
                      mu = mu, abx.s = abx.s, abx.r = abx.r, 
                      p.infect = p.infect, cum.r.1 = cum.r.1, p.r.day1 = p.r.day1, p.r.after = p.r.after,
                      short_dur = short_dur, long_dur = long_dur), 
             short_output = short_output, long_output = long_output)
  
  return(out)
}
