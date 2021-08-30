source('models/los_abx_matrix_varydur.R')
source('models/summary_los.R')
source('models/model_cocarriage5state.R')

# set transmission to 0, then low then high after 150 days each (excluding 150 days of burn in)

# get data per iteration in terms of patient matrix, los, abx and daily states 
get_output_trialdata_varytransmission_cocarriage <- function(n.bed, max.los, n.day, 
                                                             prop_R, prop_r, prop_Sr, prop_S,
                                                             bif, pi_ssr, repop.s, repop.r,
                                                             mu, abx.s, abx.r, 
                                                             p.infect, cum.r.1, p.r.day1, p.r.after,
                                                             meanDur, dur.type, timestep, iterations, burn_in){
  
  message(paste0('running ', dur.type, ' duration for ', iterations, ' iterations...'))
  
  totaldays = (n.day + 2 * (n.day - burn_in)) # set transmission to 0, then low then high after 150 days each (excluding 50 days of burn in)
  
  # set up empty matrix to record outcomes 
  outcomes.to.monitor.perday = c('prop_sR_perday', 'prop_S_perday', 'prop_s_perday', 'prop_sr_perday', 'prop_Sr_perday') # per day during los
  outcomes.to.monitor.perpatient = c('abxdur.perpatient', 'abx1dur.perpatient', 'abx2dur.perpatient')
  out.perday = array(NA, dim = c(totaldays - burn_in, iterations, length(outcomes.to.monitor.perday)),
                     dimnames = list(NULL, NULL, outcomes.to.monitor.perday))
  out.perpatient = array(NA, dim = c(iterations, length(outcomes.to.monitor.perpatient)),
                         dimnames = list(NULL, outcomes.to.monitor.perpatient))
  
  for(iter in 1:iterations){ # for each iteration 
    
    # get matrix of length of stay, abx prescribed, patients admitted
    matrixes = los_abx_table_varydur(n.bed = n.bed, n.day = totaldays, max.los = max.los, 
                                     p.infect = p.infect, p.r.day1 = p.r.day1, p.r.after = p.r.after, 
                                     cum.r.1 = cum.r.1, 
                                     meanDur = meanDur, timestep=timestep)
    patient.matrix = matrixes[[1]]
    abx.matrix = matrixes[[2]]
    los.array = summary_los(patient.matrix = patient.matrix)
    
    # starting state for all the patients admitted 
    colo.matrix = colo.table(patient.matrix = patient.matrix, los=los.array, 
                             prop_R = prop_R, prop_r=prop_r, prop_Sr=prop_Sr, prop_S=prop_S)
    
    # update values for each day 
    colo.matrix_filled_trans0_iter = nextDay(patient.matrix = patient.matrix[1:n.day,], abx.matrix = abx.matrix[1:n.day,], colo.matrix = colo.matrix[1:n.day,], 
                                             pi_ssr = pi_ssr, bif = bif, mu = mu, repop.r = repop.r, 
                                             repop.s = repop.s, abx.r=abx.r, abx.s = abx.s, timestep = timestep)
    
    colo.matrix_transLOW = colo.matrix[(n.day + 1): (n.day + n.day - burn_in),] # ensure all in the 1st row have a starting status
    colo.matrix_transLOW[1, is.na(colo.matrix_transLOW[1,])] = colo.matrix_filled_trans0_iter[nrow(colo.matrix_filled_trans0_iter),][is.na(colo.matrix_transLOW[1,])]
    colo.matrix_filled_transLOW_iter = nextDay(patient.matrix = patient.matrix[(n.day + 1): (n.day + n.day - burn_in),], 
                                               abx.matrix = abx.matrix[(n.day + 1): (n.day + n.day - burn_in),], colo.matrix = colo.matrix_transLOW, 
                                               pi_ssr = pi_ssr + 0.1, bif = bif, mu = mu, repop.r = repop.r, 
                                               repop.s = repop.s, abx.r=abx.r, abx.s = abx.s, timestep = timestep)
    
    colo.matrix_transHIGH = colo.matrix[(n.day + 1 + n.day - burn_in):totaldays, ] # ensure all in the 1st row have a starting status
    colo.matrix_transHIGH[1, is.na(colo.matrix_transHIGH[1,])] = colo.matrix_filled_transLOW_iter[nrow(colo.matrix_filled_transLOW_iter),][is.na(colo.matrix_transHIGH[1,])]
    colo.matrix_filled_transHIGH_iter = nextDay(patient.matrix = patient.matrix[(n.day + 1 + n.day - burn_in):totaldays, ], 
                                                abx.matrix = abx.matrix[(n.day + 1 + n.day - burn_in):totaldays, ], colo.matrix = colo.matrix_transHIGH, 
                                                pi_ssr = pi_ssr + 0.2, bif = bif, mu = mu, repop.r = repop.r, 
                                                repop.s = repop.s, abx.r=abx.r, abx.s = abx.s, timestep = timestep)
    colo.matrix_filled_iter = rbind(colo.matrix_filled_trans0_iter, colo.matrix_filled_transLOW_iter, colo.matrix_filled_transHIGH_iter)
    
    # remove burn-in period 
    days.after.burn.in = (burn_in + 1) : totaldays
    patient.matrix.after.burn.in = patient.matrix[days.after.burn.in,]
    pt.id.after.burn.in = unique(as.vector(patient.matrix.after.burn.in))
    abx.matrix.after.burn.in = abx.matrix[days.after.burn.in,]
    carriage.status.after.burn.in = colo.matrix_filled_iter[days.after.burn.in,]
    
    ####################
    # get summary data #
    ####################
    
    # PATIENTS 
    ### proportion of R carriers for each day of observation for this iteration
    out.perday[ , iter, 'prop_sR_perday'] = rowMeans( carriage.status.after.burn.in == 'sR')
    out.perday[ , iter, 'prop_S_perday'] = rowMeans( carriage.status.after.burn.in == 'S')
    out.perday[ , iter, 'prop_s_perday'] = rowMeans( carriage.status.after.burn.in == 'ss')
    out.perday[ , iter, 'prop_sr_perday'] = rowMeans( carriage.status.after.burn.in == 'sr')
    out.perday[ , iter, 'prop_Sr_perday'] = rowMeans( carriage.status.after.burn.in == 'Sr')
    
    ### proportion of new acquisitions each day of observation divided by the susceptibles in the previous day
    # get each patient's stay in lists 
    los.after.burn.in = rle(as.vector(patient.matrix.after.burn.in))$lengths # duration of los 
    # admission.positions.after.burn.in = c(1, cumsum(los.after.burn.in) + 1)
    # discharge.positions.after.burn.in = cumsum(los.after.burn.in)
    # # new R defined as 'sR' on the day of observation when previous day not a 'sR' AND the same patient 
    # carriage.status.per.patient.after.burn.in = split(as.vector(carriage.status.after.burn.in), rep.int(seq_along(los.after.burn.in), los.after.burn.in))
    # newR_perpatient = lapply(carriage.status.per.patient.after.burn.in, function(status.per.patient){
    #   
    #   R = c(0, which(status.per.patient == 'sR'))
    #   newR = R[(which(diff(R) > 1)) + 1]
    #   newR_binary = rep(0, length(status.per.patient))
    #   newR_binary[newR] = 1
    #   
    #   return(newR_binary)
    # })
    # newR_patient_matrix = matrix(unlist(newR_perpatient), ncol = ncol(carriage.status.after.burn.in))
    # newR_ontheday = rowSums(newR_patient_matrix == 1)[-1]
    # nonR_onthepreviousday = rowSums(newR_patient_matrix == 0) [-nrow(newR_patient_matrix)]
    # out.perday[ , iter, 'prop_newR_perday'] = c(0, newR_ontheday / nonR_onthepreviousday)
    
    # ABX USED 
    ### abx duration for each patient 
    abx.per.patient.after.burn.in = split(as.vector(abx.matrix.after.burn.in), rep.int(seq_along(los.after.burn.in), los.after.burn.in))
    abx.dur = lapply(abx.per.patient.after.burn.in, function(abx.per.patient){
      c(sum(abx.per.patient == 1), sum(abx.per.patient == 2), sum(abx.per.patient > 0))
    })
    abx.dur.df = do.call('rbind', abx.dur)
    out.perpatient[iter, 'abx1dur.perpatient'] = mean(abx.dur.df[,1])
    out.perpatient[iter, 'abx2dur.perpatient'] = mean(abx.dur.df[,2])
    out.perpatient[iter, 'abxdur.perpatient'] = mean(abx.dur.df[,3])
    
  }
  
  out = list(prop_sR_perday_perparameterset = rowMeans(out.perday[, , 'prop_sR_perday']), # average over the iterations
             prop_s_perday_perparameterset = rowMeans(out.perday[, , 'prop_s_perday']), 
             prop_Sr_perday_perparameterset = rowMeans(out.perday[, , 'prop_Sr_perday']), 
             prop_sr_perday_perparameterset = rowMeans(out.perday[, , 'prop_sr_perday']), 
             prop_S_perday_perparameterset = rowMeans(out.perday[, , 'prop_S_perday']), 
             abxdur_perpatient_perparameterset = mean(out.perpatient[,'abxdur.perpatient']),
             abx1dur_perpatient_perparameterset = mean(out.perpatient[,'abx1dur.perpatient']),
             abx2dur_perpatient_perparameterset = mean(out.perpatient[,'abx2dur.perpatient']))
  
  
  return(out)
  
}


#################### run the model using the above functions  #####################
run_trialdata_varytransmission_cocarriage5state <- function(n.bed, max.los, 
                                                            prop_R, prop_r, prop_Sr, prop_S,
                                                            bif, pi_ssr, repop.s, repop.r,
                                                            mu, abx.s, abx.r, 
                                                            p.infect, cum.r.1, p.r.day1, p.r.after,
                                                            short_dur, long_dur){
  
  message('cocarriage 5 state model initiating...')
  
  timestep = 1
  iterations = 50 # from AA tests
  n.day = 300
  burn_in = 150    # from equilibrium graphs 
  
  ############ 
  #### SHORT DURATION 
  
  short_output = get_output_trialdata_varytransmission_cocarriage(n.bed = n.bed, n.day = n.day, max.los = max.los, 
                                                                  p.infect = p.infect, p.r.day1 = p.r.day1, p.r.after = p.r.after, 
                                                                  cum.r.1 = cum.r.1, 
                                                                  meanDur = short_dur, dur.type = 'short',
                                                                  prop_R = prop_R, prop_r = prop_r, prop_Sr = prop_Sr, prop_S = prop_S, 
                                                                  pi_ssr=pi_ssr, bif=bif, mu=mu, repop.r=repop.r, 
                                                                  repop.s=repop.s, abx.r=abx.r, abx.s=abx.s, timestep=timestep, iterations = iterations, 
                                                                  burn_in = burn_in)
  
  ############ 
  #### LONG DURATION 
  long_output = get_output_trialdata_varytransmission_cocarriage(n.bed = n.bed, n.day = n.day, max.los = max.los, 
                                                                 p.infect = p.infect, p.r.day1 = p.r.day1, p.r.after = p.r.after, 
                                                                 cum.r.1 = cum.r.1, 
                                                                 meanDur = long_dur, dur.type = 'long',
                                                                 prop_R = prop_R, prop_r = prop_r, prop_Sr = prop_Sr, prop_S = prop_S, 
                                                                 pi_ssr = pi_ssr, bif = bif, mu = mu, repop.r = repop.r, 
                                                                 repop.s = repop.s, abx.r = abx.r, abx.s = abx.s, timestep = timestep, iterations = iterations, 
                                                                 burn_in = burn_in)
  
  out = list(para = c(n.bed = n.bed, max.los = max.los, 
                      prop_R = prop_R, prop_r = prop_r, prop_Sr = prop_Sr, prop_S = prop_S,
                      bif = bif, pi_ssr = pi_ssr, repop.s = repop.s, repop.r = repop.r,
                      mu = mu, abx.s = abx.s, abx.r = abx.r, 
                      p.infect = p.infect, cum.r.1 = cum.r.1, p.r.day1 = p.r.day1, p.r.after = p.r.after,
                      short_dur = short_dur, long_dur = long_dur), 
             short_output = short_output, long_output = long_output)
  
  return(out)
}
