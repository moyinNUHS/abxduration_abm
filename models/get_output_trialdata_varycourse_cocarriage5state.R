source('models/los_abx_matrix_varycourse.R')
source('models/summary_los.R')
source('models/model_cocarriage5state.R')

# get data per iteration in terms of patient matrix, los, abx and daily states 

#################### run the model using the above functions  #####################
run_trialdata_varycourse_cocarriage5state <- function(n.bed, max.los, 
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
  outcomes.to.monitor = c('totalR.continuous', 'newR.continuous', 
                          'totalR.interrupted', 'newR.interrupted')
  outcomes.to.monitor.perpatient = c('abxdur.continuous.perpatient', 'abx1dur.continuous.perpatient', 'abx2dur.continuous.perpatient', 
                                     'abxcourse.continuous.perpatient', 'abx1course.continuous.perpatient', 'abx2course.continuous.perpatient',
                                     'abxdur.interrupted.perpatient', 'abx1dur.interrupted.perpatient', 'abx2dur.interrupted.perpatient', 
                                     'abxcourse.interrupted.perpatient', 'abx1course.interrupted.perpatient', 'abx2course.interrupted.perpatient')
  out = array(NA, dim = c(iterations, length(outcomes.to.monitor)), # each row is an iter, each col is an outcome type 
              dimnames = list(NULL, outcomes.to.monitor))
  out.perpatient = array(NA, dim = c(iterations, length(outcomes.to.monitor.perpatient)), # each row is an iter, each col is an outcome type 
                         dimnames = list(NULL, outcomes.to.monitor.perpatient))
  
  for(iter in 1:iterations){ # for each iteration 
    
    # get matrix of length of stay, abx prescribed, patients admitted
    matrixes = los_abx_table_course(n.bed = n.bed, n.day = n.day, max.los=max.los, 
                                    p.infect=p.infect, abx.type = abx.type, cum.r.1=cum.r.1, 
                                    meanDur = meanDur, timestep=timestep)
    patient.matrix = matrixes[[1]]
    abx.matrix.continuous = matrixes[[2]]
    abx.matrix.interrupted = matrixes[[3]]
    los.array = summary_los(patient.matrix = patient.matrix)
    
    # starting state for all the patients admitted 
    colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                             prop_R=prop_R, prop_r=prop_r, prop_Sr=prop_Sr, prop_S=prop_S)
    
    # update values for each day 
    colo.matrix_filled_continuous_iter = nextDay(patient.matrix = patient.matrix, abx.matrix=abx.matrix.continuous, colo.matrix=colo.matrix, 
                                                 pi_ssr=pi_ssr, bif=bif, mu=mu, fitness.r=fitness.r, 
                                                 repop.s=repop.s, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
    
    colo.matrix_filled_interrupted_iter = nextDay(patient.matrix = patient.matrix, abx.matrix=abx.matrix.interrupted, colo.matrix=colo.matrix, 
                                                  pi_ssr=pi_ssr, bif=bif, mu=mu, fitness.r=fitness.r, 
                                                  repop.s=repop.s, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
    
    # remove burn-in period 
    days.after.burn.in = (burn_in + 1) : n.day
    patient.matrix.after.burn.in = patient.matrix[days.after.burn.in,]
    patient.id.after.burn.in = unique(as.vector(patient.matrix.after.burn.in))
    abx.matrix.continuous.after.burn.in = abx.matrix.continuous[days.after.burn.in,]
    carriage.status.continuous.after.burn.in = colo.matrix_filled_continuous_iter[days.after.burn.in,]
    abx.matrix.interrupted.after.burn.in = abx.matrix.interrupted[days.after.burn.in,]
    carriage.status.interrupted.after.burn.in = colo.matrix_filled_interrupted_iter[days.after.burn.in,]
    
    # get each patient's abx in lists
    los.after.burn.in = rle(as.vector(patient.matrix.after.burn.in))$lengths # duration of los 
    abx.per.patient.continuous.after.burn.in = split(as.vector(abx.matrix.continuous.after.burn.in), rep.int(seq_along(los.after.burn.in), los.after.burn.in))
    abx.per.patient.interrupted.after.burn.in = split(as.vector(abx.matrix.interrupted.after.burn.in), rep.int(seq_along(los.after.burn.in), los.after.burn.in))
    
    ####################
    # get summary data #
    ####################
    
    # PATIENTS 
    ### proportion of R carriers per iteration 
    out[iter, 'totalR.continuous'] = mean(carriage.status.continuous.after.burn.in == 'sR')
    out[iter, 'newR.continuous'] = mean(carriage.status.continuous.after.burn.in == 'sR')
    out[iter, 'totalR.interrupted'] = mean(carriage.status.interrupted.after.burn.in == 'sR')
    out[iter, 'totalR.interrupted'] = mean(carriage.status.interrupted.after.burn.in == 'sR')
    
    ### proportion of new acquisitions per iteration 
    admission.positions = (1: length(as.vector(patient.matrix.after.burn.in)))[!duplicated(as.vector(patient.matrix.after.burn.in))]
    discharge.positions = admission.positions + los.after.burn.in - 1
    out[iter, 'newR.continuous'] = sum(as.vector(carriage.status.continuous.after.burn.in)[admission.positions] != 'sR' & 
                                         as.vector(carriage.status.continuous.after.burn.in)[discharge.positions] == 'sR') / 
      sum(as.vector(carriage.status.continuous.after.burn.in)[admission.positions] != 'sR')
    
    out[iter, 'newR.interrupted'] = sum(as.vector(carriage.status.interrupted.after.burn.in)[admission.positions] != 'sR' & 
                                          as.vector(carriage.status.interrupted.after.burn.in)[discharge.positions] == 'sR') /
      sum(as.vector(carriage.status.interrupted.after.burn.in)[admission.positions] != 'sR')
    
    # ABX USED 
    ### abx duration for each patient 
    abx.dur = lapply(abx.per.patient.continuous.after.burn.in, function(abx.per.patient){
      c(sum(abx.per.patient == 1), sum(abx.per.patient == 2), sum(abx.per.patient > 0), 
        sum(rle(abx.per.patient)$values == 1), sum(rle(abx.per.patient)$values == 2), sum(rle(abx.per.patient)$values > 0))
    })
    abx.dur.df = do.call('rbind', abx.dur)
    out.perpatient[iter, 'abx1dur.continuous.perpatient'] = mean(abx.dur.df[,1], na.rm = T)
    out.perpatient[iter, 'abx2dur.continuous.perpatient'] = mean(abx.dur.df[,2], na.rm = T)
    out.perpatient[iter, 'abxdur.continuous.perpatient'] = mean(abx.dur.df[,3], na.rm = T)
    out.perpatient[iter, 'abx1course.continuous.perpatient'] = mean(abx.dur.df[,4], na.rm = T)
    out.perpatient[iter, 'abx2course.continuous.perpatient'] = mean(abx.dur.df[,5], na.rm = T)
    out.perpatient[iter, 'abxcourse.continuous.perpatient'] = mean(abx.dur.df[,6], na.rm = T)
    
    abx.dur = lapply(abx.per.patient.interrupted.after.burn.in, function(abx.per.patient){
      c(sum(abx.per.patient == 1), sum(abx.per.patient == 2), sum(abx.per.patient > 0), 
        sum(rle(abx.per.patient)$values == 1), sum(rle(abx.per.patient)$values == 2), sum(rle(abx.per.patient)$values > 0))
    })
    abx.dur.df = do.call('rbind', abx.dur)
    out.perpatient[iter, 'abx1dur.interrupted.perpatient'] = mean(abx.dur.df[,1], na.rm = T)
    out.perpatient[iter, 'abx2dur.interrupted.perpatient'] = mean(abx.dur.df[,2], na.rm = T)
    out.perpatient[iter, 'abxdur.interrupted.perpatient'] = mean(abx.dur.df[,3], na.rm = T)
    out.perpatient[iter, 'abx1course.interrupted.perpatient'] = mean(abx.dur.df[,4], na.rm = T)
    out.perpatient[iter, 'abx2course.interrupted.perpatient'] = mean(abx.dur.df[,5], na.rm = T)
    out.perpatient[iter, 'abxcourse.interrupted.perpatient'] = mean(abx.dur.df[,6], na.rm = T)
  }
  
  
  out = list(prop_R_continuous_perparameterset = mean(out[, 'totalR.continuous']),
             prop_R_interrupted_perparameterset = mean(out[, 'totalR.interrupted']),
             prop_newR_continuous_perparameterset = mean(out[, 'newR.continuous']),
             prop_newR_interrupted_perparameterset = mean(out[, 'newR.interrupted']),
             abxdur_continuous_perpatient_perparameterset = mean(out.perpatient[,'abxdur.continuous.perpatient']),
             abx1dur_continuous_perpatient_perparameterset = mean(out.perpatient[,'abx1dur.continuous.perpatient']),
             abx2dur_continuous_perpatient_perparameterset = mean(out.perpatient[,'abx2dur.continuous.perpatient']), 
             abxcourse_continuous_perpatient_perparameterset = mean(out.perpatient[,'abxcourse.continuous.perpatient']),
             abx1course_continuous_perpatient_perparameterset = mean(out.perpatient[,'abx1course.continuous.perpatient']),
             abx2course_continuous_perpatient_perparameterset = mean(out.perpatient[,'abx2course.continuous.perpatient']),
             abxdur_interrupted_perpatient_perparameterset = mean(out.perpatient[,'abxdur.interrupted.perpatient']),
             abx1dur_interrupted_perpatient_perparameterset = mean(out.perpatient[,'abx1dur.interrupted.perpatient']),
             abx2dur_interrupted_perpatient_perparameterset = mean(out.perpatient[,'abx2dur.interrupted.perpatient']), 
             abxcourse_interrupted_perpatient_perparameterset = mean(out.perpatient[,'abxcourse.interrupted.perpatient']),
             abx1course_interrupted_perpatient_perparameterset = mean(out.perpatient[,'abx1course.interrupted.perpatient']),
             abx2course_interrupted_perpatient_perparameterset = mean(out.perpatient[,'abx2course.interrupted.perpatient']))
  
  return(out)
}

