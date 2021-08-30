source('models/los_abx_matrix_varystarttime.R')
source('models/summary_los.R')
source('models/model_cocarriage5state.R')

# get data per iteration in terms of patient matrix, los, abx and daily states 

#################### run the model using the above functions  #####################
run_trialdata_varystarttime_cocarriage5state <- function(n.bed, max.los, 
                                                         prop_R, prop_r, prop_Sr, prop_S,
                                                         bif, pi_ssr, repop.s, fitness.r,
                                                         mu, abx.s, abx.r, abx.type,
                                                         p.infect, cum.r.1, meanDur){
  
  abx.type = round(abx.type)
  
  timestep = 1
  iterations = 50 # from AA tests
  n.day = 300
  burn_in = 150    # from equilibrium graphs 
  
  # set up empty matrix to record outcomes 
  outcomes.to.monitor = c('totalR.startonadmission', 'newR.startonadmission', 
                          'totalR.startlater', 'newR.startlater')
  outcomes.to.monitor.perpatient = c('abxdur.startonadmission.perpatient', 'abx1dur.startonadmission.perpatient', 'abx2dur.startonadmission.perpatient', 
                                     'abxdur.startlater.perpatient', 'abx1dur.startlater.perpatient', 'abx2dur.startlater.perpatient')
  out = array(NA, dim = c(iterations, length(outcomes.to.monitor)), # each row is an iter, each col is an outcome type 
              dimnames = list(NULL, outcomes.to.monitor))
  out.perpatient = array(NA, dim = c(iterations, length(outcomes.to.monitor.perpatient)), # each row is an iter, each col is an outcome type 
                         dimnames = list(NULL, outcomes.to.monitor.perpatient))
  
  for(iter in 1:iterations){ # for each iteration 
    
    # get matrix of length of stay, abx prescribed, patients admitted
    matrixes = los_abx_table_starttime(n.bed = n.bed, n.day = n.day, max.los=max.los, 
                                       p.infect=p.infect, abx.type = abx.type, cum.r.1=cum.r.1, 
                                       meanDur = meanDur, timestep=timestep)
    patient.matrix = matrixes[[1]]
    abx.matrix.startonadmission = matrixes[[2]]
    abx.matrix.startlater = matrixes[[3]]
    los.array = summary_los(patient.matrix = patient.matrix)
    
    # starting state for all the patients admitted 
    colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                             prop_R=prop_R, prop_r=prop_r, prop_Sr=prop_Sr, prop_S=prop_S)
    
    # update values for each day 
    colo.matrix_filled_startonadmission_iter = nextDay(patient.matrix = patient.matrix, abx.matrix=abx.matrix.startonadmission, colo.matrix=colo.matrix, 
                                                       pi_ssr=pi_ssr, bif=bif, mu=mu, fitness.r=fitness.r, 
                                                       repop.s=repop.s, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
    
    colo.matrix_filled_startlater_iter = nextDay(patient.matrix = patient.matrix, abx.matrix=abx.matrix.startlater, colo.matrix=colo.matrix, 
                                                 pi_ssr=pi_ssr, bif=bif, mu=mu, fitness.r=fitness.r, 
                                                 repop.s=repop.s, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
    
    # remove burn-in period 
    days.after.burn.in = (burn_in + 1) : n.day
    patient.matrix.after.burn.in = patient.matrix[days.after.burn.in,]
    patient.id.after.burn.in = unique(as.vector(patient.matrix.after.burn.in))
    abx.matrix.startonadmission.after.burn.in = abx.matrix.startonadmission[days.after.burn.in,]
    carriage.status.startonadmission.after.burn.in = colo.matrix_filled_startonadmission_iter[days.after.burn.in,]
    abx.matrix.startlater.after.burn.in = abx.matrix.startlater[days.after.burn.in,]
    carriage.status.startlater.after.burn.in = colo.matrix_filled_startlater_iter[days.after.burn.in,]
    
    # get each patient's abx in lists
    los.after.burn.in = rle(as.vector(patient.matrix.after.burn.in))$lengths # duration of los 
    abx.per.patient.startonadmission.after.burn.in = split(as.vector(abx.matrix.startonadmission.after.burn.in), rep.int(seq_along(los.after.burn.in), los.after.burn.in))
    abx.per.patient.startlater.after.burn.in = split(as.vector(abx.matrix.startlater.after.burn.in), rep.int(seq_along(los.after.burn.in), los.after.burn.in))
    
    ####################
    # get summary data #
    ####################
    
    # PATIENTS 
    ### proportion of R carriers per iteration 
    out[iter, 'totalR.startonadmission'] = mean(carriage.status.startonadmission.after.burn.in == 'sR')
    out[iter, 'newR.startonadmission'] = mean(carriage.status.startonadmission.after.burn.in == 'sR')
    out[iter, 'totalR.startlater'] = mean(carriage.status.startlater.after.burn.in == 'sR')
    out[iter, 'totalR.startlater'] = mean(carriage.status.startlater.after.burn.in == 'sR')
    
    ### proportion of new acquisitions per iteration 
    admission.positions = (1: length(as.vector(patient.matrix.after.burn.in)))[!duplicated(as.vector(patient.matrix.after.burn.in))]
    discharge.positions = admission.positions + los.after.burn.in - 1
    out[iter, 'newR.startonadmission'] = sum(as.vector(carriage.status.startonadmission.after.burn.in)[admission.positions] != 'sR' & 
                                               as.vector(carriage.status.startonadmission.after.burn.in)[discharge.positions] == 'sR') / 
      sum(as.vector(carriage.status.startonadmission.after.burn.in)[admission.positions] != 'sR')
    
    out[iter, 'newR.startlater'] = sum(as.vector(carriage.status.startlater.after.burn.in)[admission.positions] != 'sR' & 
                                         as.vector(carriage.status.startlater.after.burn.in)[discharge.positions] == 'sR') /
      sum(as.vector(carriage.status.startlater.after.burn.in)[admission.positions] != 'sR')
    
    # ABX USED 
    ### abx duration for each patient 
    abx.dur = lapply(abx.per.patient.startonadmission.after.burn.in, function(abx.per.patient){
      c(sum(abx.per.patient == 1), sum(abx.per.patient == 2), sum(abx.per.patient > 0))
    })
    abx.dur.df = do.call('rbind', abx.dur)
    out.perpatient[iter, 'abx1dur.startonadmission.perpatient'] = mean(abx.dur.df[,1], na.rm = T)
    out.perpatient[iter, 'abx2dur.startonadmission.perpatient'] = mean(abx.dur.df[,2], na.rm = T)
    out.perpatient[iter, 'abxdur.startonadmission.perpatient'] = mean(abx.dur.df[,3], na.rm = T)
    
    abx.dur = lapply(abx.per.patient.startlater.after.burn.in, function(abx.per.patient){
      c(sum(abx.per.patient == 1), sum(abx.per.patient == 2), sum(abx.per.patient > 0))
    })
    abx.dur.df = do.call('rbind', abx.dur)
    out.perpatient[iter, 'abx1dur.startlater.perpatient'] = mean(abx.dur.df[,1], na.rm = T)
    out.perpatient[iter, 'abx2dur.startlater.perpatient'] = mean(abx.dur.df[,2], na.rm = T)
    out.perpatient[iter, 'abxdur.startlater.perpatient'] = mean(abx.dur.df[,3], na.rm = T)
  }
  
  
  out = list(prop_R_startonadmission_perparameterset = mean(out[, 'totalR.startonadmission']),
             prop_R_startlater_perparameterset = mean(out[, 'totalR.startlater']),
             prop_newR_startonadmission_perparameterset = mean(out[, 'newR.startonadmission']),
             prop_newR_startlater_perparameterset = mean(out[, 'newR.startlater']),
             abxdur_startonadmission_perpatient_perparameterset = mean(out.perpatient[,'abxdur.startonadmission.perpatient']),
             abx1dur_startonadmission_perpatient_perparameterset = mean(out.perpatient[,'abx1dur.startonadmission.perpatient']),
             abx2dur_startonadmission_perpatient_perparameterset = mean(out.perpatient[,'abx2dur.startonadmission.perpatient']), 
             abxdur_startlater_perpatient_perparameterset = mean(out.perpatient[,'abxdur.startlater.perpatient']),
             abx1dur_startlater_perpatient_perparameterset = mean(out.perpatient[,'abx1dur.startlater.perpatient']),
             abx2dur_startlater_perpatient_perparameterset = mean(out.perpatient[,'abx2dur.startlater.perpatient']))
  
  return(out)
}

