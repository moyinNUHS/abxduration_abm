source('models/msm_util_rtnorm.R')
load("models/norm_ecdf.Rdata")

#### in this co-carriage model, the difference between the two wards is 1 continuous course vs many short courses 
### such that 
# able to define which antibiotic is given 
# vary continuous vs interrupted antibiotics


###################
# define los and antibiotic duration  #####################
los.abx.table.course <- function(n.bed, n.day, max.los, abx.type, 
                                 p.infect, cum.r.1,
                                 meanDur, timestep){
  
  abx.type = round(abx.type)
  
  # maximum number of patients that can be admitted to the ward in n.day
  # (maximum number of patients if everyone stays only for 1 day)
  n.patient = ceiling(n.bed * n.day)
  patient.id = 1:n.patient
  
  # generate length of stay for each patient IF no hospital acquired infection 
  # (truncated normal distribution)
  all_los = ceiling(rexp(n.patient, 1/max.los))
  all_los[all_los > (5 * max.los)] = ceiling(max.los)
  
  ####################
  # DETERMINE ABX COURSES 
  
  # first decide if patients are on antibiotics on day 1 of admission 
  #number of broad vs narrow spectrum antibiotics on day 1 of admission 
  no.abx.day1 = round(n.patient * p.infect) #number of antibiotics prescribed for community acquired infection 
  #identify patient ids who were prescribed antibiotics on day 1 of admission
  id.abx.day1 = sample(patient.id, no.abx.day1) #id of patients prescribed antibiotics on day 1
  
  # then decide if going to receive antibiotics after admission - risk of acquiring HAI increases with length of stay 
  # (cumulative normal distribution)
  # if prescribed antibiotics after admission then prolong admission to complete course 
  prob.r.after = norm_ecdf((1/cum.r.1) * (1 : max(all_los))) # get the cumulative probability for each day
  # norm_ecdf - imported normalized function
  abx.dayafter = lapply(patient.id, function(id){
    # if the patient will receive antibiotics after admission - 1 if will be prescribed antibiotics on the day during los 
    abx.dayafter.binary = rbinom(all_los[id], 1, prob = prob.r.after[1:all_los[id]])
    # start dates if given antibiotics
    if (sum(abx.dayafter.binary) > 0) {which(abx.dayafter.binary == 1)} else { 0 }
  })
  id.abx.dayafter = patient.id[which(lengths(abx.dayafter) > 1)]
  
  #allocate antibiotics for the patients since admission
  all_abx = lapply(patient.id, function(id){ # for each patient 
    
    # create empty vector with 0 with length same as length of stay
    fill = rep(0, all_los[id])
    
    # fill antibiotics on admission
    if (length(id.abx.day1) > 0 & id %in% id.abx.day1) { # if the patient is given antibiotics on admission
      # antibiotic duration randomly drawn from a truncated normal distribution
      dur = round(rtnorm(1, mean = meanDur, lower = 1))
      fill[1: dur] = abx.type
    } 
    
    # fill antibiotics after admission
    if (length(id.abx.dayafter) > 0 & id %in% id.abx.dayafter) {
      
      courses = abx.dayafter[[id]][which(abx.dayafter[[id]] > 0)]
      no.courses = length(courses)
      dur = round(rtnorm(no.courses, mean = meanDur, lower = 1))
      
      for (i in 1:no.courses){
        start = courses[i]
        fill[start : (start + dur[i] - 1)] = abx.type
      }
      
    }
    
    # get continuous vs courses
    if (sum(fill) > 0) {
      start = min(which(fill == abx.type))
      total.dur = sum(fill == abx.type)
      fill.continuous = rep(0, length(fill))
      fill.continuous[start : (start + total.dur - 1)] = abx.type
      
      fill.interrupted = sample(fill)
    } else {
      fill.continuous = fill.interrupted = fill
    }
    
    return(list(fill.continuous = fill.continuous, fill.interrupted = fill.interrupted))
  })
  
  all_abx_continuous = lapply(all_abx, `[[`, 1)
  all_abx_interrupted = lapply(all_abx, `[[`, 2)
  
  ####################
  #PUT ABX COURSES AND LOS INTO PATIENT and ABX MATRIX
  #make up a matrix of number of days we want to observe (rows) -
  #against number of beds in the ward (columns)
  all_los = lengths(all_abx_continuous) #update length of stay incorporating abx given after admission
  patient.matrix = abx.matrix.continuous = abx.matrix.interrupted = matrix(NA, nrow = n.day, ncol = n.bed)
  idx = 1
  for(j in 1:n.bed){
    
    sum_los = cumsum(all_los[idx:length(all_los)]) #cumulative stay of patients 
    
    #find last patient in the last row (last day of observation)
    end.pt.id = suppressWarnings(max(which(sum_los < n.day)) + idx)
    if (end.pt.id == -Inf) {end.pt.id = idx}
    if (end.pt.id > length(all_los)) {end.pt.id = length(all_los)}
    #fill matrix
    patient.matrix[, j] = rep(idx:end.pt.id, all_los[idx: end.pt.id])[1:n.day]
    abx.matrix.continuous[,j] = unlist(all_abx_continuous[idx:end.pt.id])[1:n.day]
    abx.matrix.interrupted[,j] = unlist(all_abx_interrupted[idx:end.pt.id])[1:n.day]
    
    idx = end.pt.id + 1
  }
  
  patient.matrix = patient.matrix[rep(1:n.day, each = timestep), ]
  abx.matrix.continuous = abx.matrix.continuous[rep(1:n.day, each = timestep), ]
  abx.matrix.interrupted = abx.matrix.interrupted[rep(1:n.day, each = timestep), ]
  
  return(list(patient.matrix, abx.matrix.continuous, abx.matrix.interrupted))
}

#frequency summary of patient.matrix - patient id against number of days of stay for each patient
summary.los <- function(patient.matrix){  
  
  # Summarize how often each patient.id is encountered to get days for each id
  los.dur = table(patient.matrix) ##SLOW
  los_duration = array(dim = c(2, max(patient.matrix)))
  # Attach patient ID on 1st row
  los_duration[1,] = 1:max(patient.matrix)
  # Put summary of days on 2nd row
  los_duration[2,] = los.dur
  
  return(los_duration)
}

###################
# define initial states of each patient  #####################
colo.table <- function(patient.matrix, los.array, prop_R, prop_r, prop_Sr, prop_S){
  
  prob_start_sR = prop_R
  #S, ss, Sr, sr come from (1-prop_R)
  #Sr, sr come from prop_r*(1-prop_R)
  prob_start_Sr = prop_Sr * prop_r * (1 - prop_R)
  prob_start_sr = (1 - prop_Sr) * prop_r * (1 - prop_R)
  #S, ss come from (1-prop_r)*(1-prop_R)
  prob_start_S = prop_S * (1 - prop_r) * (1 - prop_R)
  prob_start_ss = (1 - prop_S) * (1 - prop_r) * (1 - prop_R)
  
  prob_StartBact_bi = c(prob_start_S, prob_start_Sr, prob_start_sR, prob_start_sr)
  
  #Generating a vector of random status with runif (change for other distribution)
  number_of_patients = dim(los.array)[2]
  Patient_unif = runif(number_of_patients,0,1)
  Patient_StartBact = rep(NA, number_of_patients)
  Patient_StartBact[Patient_unif > sum(prob_StartBact_bi)] = 'ss'
  Patient_StartBact[(Patient_unif > sum(prob_StartBact_bi[1:3])) & (Patient_unif <= sum(prob_StartBact_bi))] = 'sr'
  Patient_StartBact[(Patient_unif > sum(prob_StartBact_bi[1:2])) & (Patient_unif <= sum(prob_StartBact_bi[1:3]))] = 'sR'
  Patient_StartBact[(Patient_unif > sum(prob_StartBact_bi[1])) & (Patient_unif <= sum(prob_StartBact_bi[1:2]))] = 'Sr'
  Patient_StartBact[Patient_unif <= prob_StartBact_bi[1]] = 'S'
  
  #Creating array for carriage status
  array_StartBact = matrix(NA, nrow = nrow(patient.matrix), ncol = ncol(patient.matrix))
  
  # Fill generated bacterial in the first day of each patient entering the ward
  end_idx = 1
  for(i in 1:number_of_patients){
    array_StartBact[end_idx:(end_idx + los.array[2, i] - 1)] = c(Patient_StartBact[i], rep(NA, los.array[2, i]-1))
    end_idx = end_idx + los.array[2, i]
  }
  
  return(array_StartBact)
}

#################### Update values for every day  #####################
nextDay <- function(patient.matrix, abx.matrix, colo.matrix, 
                    pi_ssr, bif, mu, repop.r, repop.s, abx.r, abx.s, timestep){
  
  #number of beds
  n.bed = ncol(patient.matrix)
  
  # adjust probabilities based on timestep
  if (timestep > 1) {
    pi_ssr = 1 - (1 - pi_ssr) ^ (1 / timestep)
    mu = 1 - (1 - mu) ^ (1 / timestep)
    repop.r = 1 - (1 - repop.r) ^ (1 / timestep)
    repop.s = 1 - (1 - repop.s) ^ (1 / timestep)
    abx.r = 1 - (1 - abx.r) ^ (1 / timestep)
    abx.s = 1 - (1 - abx.s) ^ (1 / timestep)
  }
  
  pi_Sr = pi_ssr * (1 - bif)   # pi_ssr = probability of R transmitting to S to become Sr 
  
  # For each day
  for(i in 2:nrow(patient.matrix)){ # start from 2 because for each day (first day should be filled)
    # Get the previous row subset of the entire matrix which represents previous day or timestep
    prev_step = colo.matrix[i-1, ]
    # Get the indices which are already filled by a patient entering the ward for exclusion
    already_filled = which(!is.na(colo.matrix[i, ]))
    # Get the previous row subset of the entire abx matrix which represents previous day or timestep
    prev_abx = abx.matrix[i-1, ]
    # count if there are any sR on the previous day
    R.previousday = which(prev_step == "sR")
    r_num = length(R.previousday)
    
    # probability of being transmitted R
    infect_Sr = 1 - ((1 - pi_Sr) ^ (r_num / n.bed))
    infect_ssr = 1 - ((1 - pi_ssr) ^ (r_num / n.bed))
    
    if(repop.s + infect_ssr > 1){
      stop(paste("Error: repop.s + infect_ssr > 1 in those with no antibiotics on previous day"))
    }
    
    if(repop.s + mu > 1){
      stop(paste("Error: repop.s + mu > 1 in those with no antibiotics on previous day"))
    }
    
    if(abx.s + infect_Sr > 1){
      stop(paste("Error: abx.s + infect_Sr >1 in those with narrow antibiotics on previous day"))
    }
    
    if(abx.s + infect_Sr > 1){
      stop(paste("Error: abx.s + infect_Sr > 1 in those with broad antibiotics on previous day"))
    }
    
    # First scenario: those with no antibiotics on previous day 
    ##########################################
    id_noabx = which(prev_abx == 0)
    
    # Update S
    # Get the column indices which contain S in the previous day
    S = id_noabx[id_noabx %in% which(prev_step == "S")]
    # Remove column indices that already have a starting bacterial state filled in
    S = S[!(S %in% already_filled)]
    # if there is any S (number of S > 0) in the previous day
    if(length(S)){
      # Roll a random number for each S on the previous day for clearance
      roll= runif(length(S), 0, 1)
      
      Sr_idx = S[roll < infect_Sr]
      same_idx= S[roll >= infect_Sr]
      
      colo.matrix[i, Sr_idx] = "Sr"
      colo.matrix[i, same_idx] = "S"
    }
    
    # Update ss
    ss = id_noabx[id_noabx %in% which(prev_step == "ss")]
    ss = ss[!(ss %in% already_filled)]
    if(length(ss)){
      
      roll = runif(length(ss), 0, 1)
      
      ssr_idx = ss[roll < infect_ssr]
      s_idx = ss[roll >= infect_ssr & roll < repop.s + infect_ssr]
      same_idx = ss[!(ss %in% c(ssr_idx, s_idx))]
      
      colo.matrix[i, s_idx] = "S"
      colo.matrix[i, ssr_idx] = "sr"
      colo.matrix[i, same_idx] = "ss"
    }
    
    # Update Sr
    Sr = id_noabx[id_noabx %in% which(prev_step == "Sr")]
    Sr = Sr[!(Sr %in% already_filled)]
    if(length(Sr)){
      
      # Roll a random number for each S on the previous day for clearance
      roll= runif(length(Sr), 0, 1)
      
      s_idx = Sr[roll < mu]
      same_idx= Sr[roll >= mu]
      
      colo.matrix[i, s_idx] = "S"
      colo.matrix[i, same_idx] = "Sr"
    }
    
    # Update sr
    sr = id_noabx[id_noabx %in% which(prev_step == "sr")]
    sr = sr[!(sr %in% already_filled)]
    if(length(sr)){
      
      roll= runif(length(sr), 0, 1)
      
      # if roll does not pass repop.r event and passes repop.s 
      Sr_idx = sr[roll < repop.s]
      # if roll does not pass repop.r and repop.s, and passes mu
      ss_idx = sr[roll >= repop.s & roll < repop.s + mu]
      
      same_idx= sr[!(sr %in% c(Sr_idx, ss_idx))]
      
      colo.matrix[i, Sr_idx] = "Sr"
      colo.matrix[i, ss_idx] = "ss"
      colo.matrix[i, same_idx] = "sr"
    }
    
    # Update sR
    sR = id_noabx[id_noabx %in% which(prev_step == "sR")]
    sR = sR[!(sR %in% already_filled)]
    if(length(sR)){
      
      roll = runif(length(sR), 0, 1)
      
      ssr_idx = sR[roll < mu]
      same_idx = sR[roll >= mu]
      
      colo.matrix[i, ssr_idx] = "sr"
      colo.matrix[i, same_idx] = "sR"
    }
    
    # Second scenario: those with antibiotics.s on previous day 
    ##########################################
    id_narrowabx= which(prev_abx==1)
    
    # Update S
    # Get the column indices which contain S in the previous day
    S = id_narrowabx[id_narrowabx %in% which(prev_step == "S")]
    # Remove column indices that already have a starting bacterial state filled in
    S = S[!(S %in% already_filled)]
    # if there is any S (number of S > 0) in the previous day
    if(length(S)){
      
      # Roll a random number for each S on the previous day for clearance
      roll= runif(length(S), 0, 1)
      
      Sr_idx = S[roll < infect_Sr]
      ss_idx = S[roll >= infect_Sr & roll < abx.s + infect_Sr]
      same_idx= S[!(S %in% c(Sr_idx, ss_idx))]
      
      colo.matrix[i, ss_idx] = "ss"
      colo.matrix[i, Sr_idx] = "Sr"
      colo.matrix[i, same_idx] = "S"
    }
    
    # Update ss
    ss = id_narrowabx[id_narrowabx %in% which(prev_step == "ss")]
    ss = ss[!(ss %in% already_filled)]
    if(length(ss)){
      
      roll = runif(length(ss), 0, 1)
      
      ssr_idx = ss[roll < infect_ssr]
      same_idx = ss[roll >= infect_ssr]
      
      colo.matrix[i, ssr_idx] = "sr"
      colo.matrix[i, same_idx] = "ss"
    }
    
    # Update Sr
    Sr = id_narrowabx[id_narrowabx %in% which(prev_step == "Sr")]
    Sr = Sr[!(Sr %in% already_filled)]
    if(length(Sr)){
      # Roll a random number for each Sr 
      roll= runif(length(Sr), 0, 1)
      
      ssr_idx = Sr[roll < abx.s]
      same_idx= Sr[roll >= abx.s]
      
      colo.matrix[i, ssr_idx] = "sr"
      colo.matrix[i, same_idx] = "Sr"
    }
    
    # Update sr
    ssr = id_narrowabx[id_narrowabx %in% which(prev_step == "sr")]
    ssr = ssr[!(ssr %in% already_filled)]
    if(length(ssr)){
      roll= runif(length(ssr), 0, 1)
      
      sR_idx = ssr[roll < repop.r]
      same_idx= ssr[roll >= repop.r]
      
      colo.matrix[i, sR_idx] = "sR"
      colo.matrix[i, same_idx] = "sr"
    }
    
    # Update sR 
    sR = id_narrowabx[id_narrowabx %in% which(prev_step == "sR")]
    sR = sR[!(sR %in% already_filled)]
    if(length(sR)){
      colo.matrix[i, sR] = "sR"
    }
    
    # Third scenario: those with antibiotics.r on previous day 
    ##########################################
    id_broadabx= which(prev_abx==2)
    
    # Update S
    # Get the column indices which contain S in the previous day
    S = id_broadabx[id_broadabx %in% which(prev_step == "S")]
    # Remove column indices that already have a starting bacterial state filled in
    S = S[!(S %in% already_filled)]
    # if there is any S (number of S > 0) in the previous day
    if(length(S)){
      # Roll a random number for each S on the previous day for clearance
      roll= runif(length(S), 0, 1)
      
      Sr_idx = S[roll < infect_Sr]
      ss_idx = S[roll >= infect_Sr & roll < abx.s + infect_Sr]
      same_idx= S[!(S %in% c(Sr_idx, ss_idx))]
      
      colo.matrix[i, ss_idx] = "ss"
      colo.matrix[i, Sr_idx] = "Sr"
      colo.matrix[i, same_idx] = "S"
    }
    
    # Update ss
    ss = id_broadabx[id_broadabx %in% which(prev_step == "ss")]
    ss = ss[!(ss %in% already_filled)]
    if(length(ss)){
      roll = runif(length(ss), 0, 1)
      
      ssr_idx = ss[roll < infect_ssr]
      same_idx = ss[roll >= infect_ssr]
      
      colo.matrix[i, ssr_idx] = "sr"
      colo.matrix[i, same_idx] = "ss"
    }
    
    # Update Sr
    Sr = id_broadabx[id_broadabx %in% which(prev_step == "Sr")]
    Sr = Sr[!(Sr %in% already_filled)]
    if(length(Sr)){
      # Roll a random number for each S on the previous day for clearance
      roll= runif(length(Sr), 0, 1)
      
      ssr_idx = Sr[roll < abx.s]
      same_idx= Sr[roll >= abx.s]
      
      colo.matrix[i, ssr_idx] = "sr"
      colo.matrix[i, same_idx] = "Sr"
    }
    
    # Update sr
    sr = id_broadabx[id_broadabx %in% which(prev_step == "sr")]
    sr = sr[!(sr %in% already_filled)]
    if(length(sr)){
      roll= runif(length(sr), 0, 1)
      
      ss_idx = sr[roll < abx.r]
      sR_idx = sr[roll >= abx.r & roll< abx.r + repop.r]
      same_idx= sr[roll >= abx.r + repop.r]
      
      colo.matrix[i, ss_idx] = "ss"
      colo.matrix[i, sR_idx] = "sR"
      colo.matrix[i, same_idx] = "sr"
    }
    
    # Update sR
    sR = id_broadabx[id_broadabx %in% which(prev_step == "sR")]
    sR = sR[!(sR %in% already_filled)]
    if(length(sR)){
      roll = runif(length(sR), 0, 1)
      
      ssr_idx = sR[roll < abx.r]
      same_idx= sR[roll >= abx.r]
      
      colo.matrix[i, ssr_idx] = "sr"
      colo.matrix[i, same_idx] = "sR"
    }
  }
  
  return(colo.matrix)
}



#################### run the model using the above functions  #####################

run_model_cocarriage5state_course <- function(n.bed, max.los, 
                                              prop_R, prop_r, prop_Sr, prop_S,
                                              bif, pi_ssr, repop.s, repop.r,
                                              mu, abx.s, abx.r, abx.type,
                                              p.infect, cum.r.1, meanDur){
  
  abx.type = round(abx.type)
  
  message('cocarriage 5 state model initiating...')
  
  timestep = 1
  iterations = 50 # from AA tests
  n.day = 300
  burn_in = 150    # from equilibrium graphs 
  
  # empty matrix to store output - each row is a day, each col is a iteration 
  iter_totalsR_continuous = iter_totalsR_interrupted = matrix(NA, nrow = n.day, ncol = iterations)
  iter_newsR_continuous = iter_newsR_interrupted = c()
  iter_abx_continuous_courses = iter_abx_interrupted_courses = c()
  iter_abx_continuous_days = iter_abx_interrupted_days = c()
  iter_abxdays =  c()
  iter_abxcourses = c()
  
  for(iter in 1:iterations){ # for each iteration 
    
    # get matrix of length of stay, abx prescribed, patients admitted
    matrixes = los.abx.table.course(n.bed = n.bed, n.day = n.day, max.los=max.los, 
                                    p.infect=p.infect, abx.type = abx.type, cum.r.1=cum.r.1, 
                                    meanDur = meanDur, timestep=timestep)
    patient.matrix = matrixes[[1]]
    abx.matrix.continuous = matrixes[[2]]
    abx.matrix.interrupted = matrixes[[3]]
    los.array = summary.los(patient.matrix = patient.matrix)
    
    # starting state for all the patients admitted 
    colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                             prop_R=prop_R, prop_r=prop_r, prop_Sr=prop_Sr, prop_S=prop_S)
    
    # update values for each day 
    colo.matrix_filled_continuous_iter = nextDay(patient.matrix = patient.matrix, abx.matrix=abx.matrix.continuous, colo.matrix=colo.matrix, 
                                                 pi_ssr=pi_ssr, bif=bif, mu=mu, repop.r=repop.r, 
                                                 repop.s=repop.s, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
    
    colo.matrix_filled_interrupted_iter = nextDay(patient.matrix = patient.matrix, abx.matrix=abx.matrix.interrupted, colo.matrix=colo.matrix, 
                                                  pi_ssr=pi_ssr, bif=bif, mu=mu, repop.r=repop.r, 
                                                  repop.s=repop.s, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
    
    # summary of the output 
    ### sRs per day 
    total_sR_perday = matrix(rowSums(colo.matrix_filled_continuous_iter == "sR"), ncol = timestep, byrow=T)
    iter_totalsR_continuous[, iter]= rowMeans(total_sR_perday)
    
    total_sR_perday = matrix(rowSums(colo.matrix_filled_interrupted_iter == "sR"), ncol = timestep, byrow=T)
    iter_totalsR_interrupted[, iter]= rowMeans(total_sR_perday)
    
    ### New R at discharge when not admitted as R 
    patient.matrix.exlude.burnin = patient.matrix[((burn_in + 1) : n.day), ]   # patient matrix after the first `burn_in` days 
    patient.id.exlude.burnin = unique(as.vector(patient.matrix.exlude.burnin)) # patient ids admitted after the first `burn_in` days 
    admit.positions = c(1, cumsum(los.array[2,]) + 1)[1:length(los.array[2,])]
    admit.positions.exlude.burnin = admit.positions[patient.id.exlude.burnin]  # admission days of patients admitted after the first `burn_in` days 
    discharge.positions = cumsum(los.array[2,])
    discharge.positions.exlude.burnin = discharge.positions[patient.id.exlude.burnin] # discharge days of patients admitted after the first `burn_in` days 
    
    new_sR_peradmission = lapply(list(colo.matrix_filled_continuous_iter, colo.matrix_filled_interrupted_iter), function(x){
      colo_vector_filled_iter = as.vector(x)
      if ('sR' %in%  x) {
        new_sR = length(which(colo_vector_filled_iter[admit.positions.exlude.burnin] != 'sR' & colo_vector_filled_iter[discharge.positions.exlude.burnin] == 'sR'))
        new_sR_peradmission = new_sR/ length(admit.positions.exlude.burnin)
      } else {
        new_sR_peradmission = 0
      }
      return(new_sR_peradmission)
    })
    
    iter_newsR_continuous[iter] = new_sR_peradmission[[1]]
    iter_newsR_interrupted[iter] = new_sR_peradmission[[2]]
    
    #keep record of days and number of courses of abx prescribed for each patient 
    abx.matrix.continuous.exlude.burnin = abx.matrix.continuous[((burn_in + 1) : n.day), ] # abx matrix after the first `burn_in` days 
    abx.matrix.interrupted.exlude.burnin = abx.matrix.interrupted[((burn_in + 1) : n.day), ] # abx matrix after the first `burn_in` days 
    
    los.exclude.burnin = rle(as.vector(patient.matrix.exlude.burnin))$lengths
    
    abx.continuous.per.patient.exclude.burnin = split(as.vector(abx.matrix.continuous.exlude.burnin), rep.int(seq_along(los.exclude.burnin), los.exclude.burnin))
    abx.interrupted.per.patient.exclude.burnin = split(as.vector(abx.matrix.interrupted.exlude.burnin), rep.int(seq_along(los.exclude.burnin), los.exclude.burnin))
    
    abx.continuous.per.patient.exclude.burnin.counts = lapply(abx.continuous.per.patient.exclude.burnin, function(x){ # for each patient
      
      counts = rle(x)
      
      abx.days = sum(counts$lengths[which(counts$values == abx.type)])
      if (length(abx.days) == 0) {abx.days = 0}
      abx.courses = sum(counts$values == abx.type)
      
      return(c(abx.continuous.days = abx.days, abx.continuous.courses = abx.courses))
    })
    
    abx.interrupted.per.patient.exclude.burnin.counts = lapply(abx.interrupted.per.patient.exclude.burnin, function(x){ # for each patient
      
      counts = rle(x)
      
      abx.days = sum(counts$lengths[which(counts$values == abx.type)])
      if (length(abx.days) == 0) {abx.days = 0}
      abx.courses = sum(counts$values == abx.type)
      
      return(c(abx.interrupted.days = abx.days, abx.interrupted.courses = abx.courses))
    })
    
    abx.per.patient.exclude.burnin.counts.df = cbind(do.call('rbind', abx.continuous.per.patient.exclude.burnin.counts), 
                                                     do.call('rbind', abx.interrupted.per.patient.exclude.burnin.counts))
    
    iter_abx_continuous_days[iter] = mean(abx.per.patient.exclude.burnin.counts.df[,'abx.continuous.days']) # mean per patient 
    iter_abx_interrupted_days[iter] = mean(abx.per.patient.exclude.burnin.counts.df[,'abx.interrupted.days'])
    iter_abx_continuous_courses[iter] = mean(abx.per.patient.exclude.burnin.counts.df[,'abx.continuous.courses'])
    iter_abx_interrupted_courses[iter] = mean(abx.per.patient.exclude.burnin.counts.df[,'abx.interrupted.courses'])
  }
  
  # Discard first `burn_in` days as burn-in
  total_sR_periter_perday = rowSums(iter_totalsR_continuous[(burn_in + 1) : nrow(iter_totalsR_continuous), , drop = FALSE]) / iterations / n.bed
  totalsR_continuous = mean(total_sR_periter_perday) # mean of the days 
  
  total_sR_periter_perday = rowSums(iter_totalsR_interrupted[(burn_in + 1) : nrow(iter_totalsR_interrupted), , drop = FALSE]) / iterations / n.bed
  totalsR_interrupted = mean(total_sR_periter_perday) # mean of the days 
  
  newsR_continuous = mean(iter_newsR_continuous) # mean per iteration
  newsR_interrupted = mean(iter_newsR_interrupted) # mean per iteration
  
  out = c(totalR_continuous = totalsR_continuous, totalR_interrupted = totalsR_interrupted, 
          newR_continuous = newsR_continuous, newR_interrupted = newsR_interrupted, 
          abx_continuous_days = mean(iter_abx_continuous_days), abx_interrupted_days = mean(iter_abx_interrupted_days), 
          abx_continuous_courses = mean(iter_abx_continuous_courses), abx_interrupted_courses = mean(iter_abx_interrupted_courses))
  
  print(out)
  
  return(out)
}

