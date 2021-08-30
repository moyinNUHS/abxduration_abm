source('~/Documents/nBox/git_projects/indiv_abxduration_incmetaanalysis/indiv_abxduration/models/msm_util_rtnorm.R')

## Co-carriage 5 state model 
# allows co-carriage of R and S bacteria 

#################### define initial states of each patient  #####################
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
                    pi_ssr, bif, mu, repop.s, fitness.r,
                    abx.r, abx.s, timestep){
  
  #number of beds
  n.bed = ncol(patient.matrix)
  
  # rate of R growth 
  pop.r = repop.s * fitness.r 
  
  # adjust probabilities based on timestep
  if (timestep > 1) {
    pi_ssr = 1 - (1 - pi_ssr) ^ (1 / timestep)
    mu = 1 - (1 - mu) ^ (1 / timestep)
    pop.r = 1 - (1 - pop.r) ^ (1 / timestep)
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
      
      roll = runif(length(sr), 0, 1)
      
      # if roll passes repop.s 
      Sr_idx = sr[roll < repop.s]
      # if roll does not pass repop.s, and passes mu
      ss_idx = sr[roll >= repop.s & roll < repop.s + mu]
      # if roll does not pass repop.s, mu and passes pop.r
      sR_idx = sr[roll >= repop.s + mu & roll < repop.s + mu + pop.r]
      
      same_idx= sr[!(sr %in% c(Sr_idx, ss_idx, sR_idx))]
      
      colo.matrix[i, Sr_idx] = "Sr"
      colo.matrix[i, ss_idx] = "ss"
      colo.matrix[i, sR_idx] = "sR"
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
    
    # Second scenario: those with antibiotic.s on previous day 
    ##########################################
    id_narrowabx = which(prev_abx == 1)
    
    # Update S
    # Get the column indices which contain S in the previous day
    S = id_narrowabx[id_narrowabx %in% which(prev_step == "S")]
    # Remove column indices that already have a starting bacterial state filled in
    S = S[!(S %in% already_filled)]
    # if there is any S (number of S > 0) in the previous day
    if(length(S)){
      
      # Roll a random number for each S on the previous day for clearance
      roll = runif(length(S), 0, 1)
      
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
      
      sR_idx = ssr[roll < pop.r]
      same_idx= ssr[roll >= pop.r]
      
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
    id_broadabx = which(prev_abx==2)
    
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
      sR_idx = sr[roll >= abx.r & roll< abx.r + pop.r]
      same_idx= sr[roll >= abx.r + pop.r]
      
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


