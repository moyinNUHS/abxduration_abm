source('~/Documents/nBox/git_projects/indiv_abxduration_incmetaanalysis/indiv_abxduration/models/msm_util_rtnorm.R')

## Population growth model 
# R and S bacteria explicitly calculated at each time step 


#distribution of R bacteria in the carriers
r_beta = rbeta(5000, 0.2, 2) 
r_beta_norm = (r_beta-min(r_beta)) / (max(r_beta) - min(r_beta)) #normalised from 0 to 1
# shape 1 and shape 2  based on rene gut data 
# d = read.csv('gutdata/Ini_CTXm_copies_qPCR.csv')
# hist(d$ini_CTXm_copies)
# hist(rbeta(10000, 0.1, 2))
# hist(r_beta_norm*exp(16))


#################### define initial states of each patient  #####################
colo.table <- function(patient.matrix, los.array, total_prop, prop_R, r_thres, K){
  
  n.day = nrow(patient.matrix)
  n.bed = ncol(patient.matrix)
  
  number_of_patients = dim(los.array)[2]
  
  #capacity for enterobacteriaceae growth (log)
  #log of the capacity is normal in distribution from rene 
  #d=read.csv('gutdata/Ini_CTXm_copies_qPCR.csv')
  #hist(d$ini_16S_log)
  total_capacity = rnorm(number_of_patients, mean = log(K)) #in log 
  total_capacity_matrix = matrix(rep(total_capacity, los.array[2,]), byrow = F, ncol = ncol(patient.matrix)) #in log 
  
  #existing population 
  #existing population mean is proportion (total_prop) of total capacity instead of a proportion of the 
  #distribution of the total capacity so that some starts at full capacity while others are not 
  #which is similar to model 2 where there is sr 
  total_existing_mean = log(total_prop * exp(total_capacity)) #in log 
  total_existing = rnorm(number_of_patients, mean = total_existing_mean) #total number of Enterobacteriaceae is a proportion of the capacity (log)
  total_existing [which(total_existing >= total_capacity)] = total_capacity [which(total_existing >= total_capacity)] #those exceeding total capacity will be given their own full capacity
  
  #amount of S and R carried 
  r_thres_log = log(r_thres * exp(total_existing)) # a proportion of the total gut capacity 
  r_thres_matrix = matrix(rep(r_thres_log, los.array[2,]), byrow = F, ncol = ncol(patient.matrix)) #in log 
  
  morethan_rthres = which(total_existing > r_thres_log)  #pick patients who have existing Enterobacteriaceae more than r_thres (log) 
  r.id = sample(morethan_rthres, size = floor(prop_R * number_of_patients)) 
  if (length(r.id) > 0) {s.id = (1:number_of_patients) [-r.id]} else {s.id = (1:number_of_patients)}
  
  r_bact = rep(NA, number_of_patients)
  # R individuals carry R 
  for (ind in r.id){  #in log 
    #expand r_beta to each individual's existing population
    r_list_full = r_beta_norm * total_existing [ind] 
    r_bact[ind] = sample(r_list_full[which(r_list_full >= r_thres_log[ind])], 1)# list of r amounts to sample from for R
  }
  
  # S individuals also carry some R (similar to Sr, sr, sR in model 2)
  for (ind in s.id){  #in log 
    #expand r_beta to each individual's existing population
    r_list_full = r_beta_norm * total_existing [ind]
    r_bact[ind] = sample(r_list_full[which(r_list_full < r_thres_log[ind])], 1)#list of r amounts to sample from for S
  }
  
  r_bact[which(r_bact<0)] = -Inf #those with less than 1 bacteria round to 0 (log)
  
  #check below equal to prop_R
  #sum(r_bact>r_thres)/number_of_patients
  #check shape 
  #hist(exp(r_bact))
  s_bact = exp(total_existing) - exp(r_bact) #abs
  s_bat_abs = ifelse(s_bact<0, 0, s_bact)
  s_bact = log(s_bat_abs) #total number of S for each patient (log)
  
  S_Bactlevelstart = matrix(NA, n.day, n.bed)
  R_Bactlevelstart = matrix(NA, n.day, n.bed)
  
  # pad with NAs
  end_idx = 1
  for(i in 1:number_of_patients){
    S_Bactlevelstart[end_idx:(end_idx + los.array[2, i] - 1)] = c(s_bact[i], rep(NA, los.array[2, i]-1))
    R_Bactlevelstart[end_idx:(end_idx + los.array[2, i] - 1)] = c(r_bact[i], rep(NA, los.array[2, i]-1))
    end_idx = end_idx + los.array[2, i]
  }
  
  return(list(S_Bactlevelstart, R_Bactlevelstart, total_capacity_matrix, r_thres_matrix)) # in log
}

#################### Update values for every day  #####################
nextDay <- function(patient.matrix, los.array, abx.matrix, colo.matrix, 
                    pi_ssr, total_prop, K, fitness.r, r_thres, r_trans, s_growth,
                    abx.s, abx.r, timestep){
  
  if (timestep > 1) {
    pi_ssr = 1 - (1 - pi_ssr) ^ (1 / timestep)
  }
  
  n.bed = ncol(patient.matrix)
  
  S_table = exp(colo.matrix[[1]]) #in abs
  R_table = exp(colo.matrix[[2]]) #in abs
  
  #capacity matrix for enterobacteriaceae growth (abs)
  total_capacity = exp(colo.matrix[[3]]) 
  
  # r_thres (abs)
  r_thres_matrix = exp(colo.matrix[[4]])
  
  # r_growth 
  r_growth = s_growth * fitness.r
  
  # For each day (first day should be filled)
  for(i in 2:nrow(patient.matrix)){
    
    # calculate how many people has R above r_thres (log) in the previous day
    r_num = sum(R_table[i-1,] >= r_thres_matrix[i-1,])
    
    if (is.na(r_num)) {r_num = 0}
    
    # from number of r, calculate probability of transmission
    prop_r = 1 - ((1 - pi_ssr) ^ (r_num / n.bed)) 
    
    ###### Convert all log scale parameters into normal scale for addition, then convert back to log
    #for each person:
    for(j in 1:ncol(patient.matrix)){
      
      #print(c(i,j))
      
      if(is.na(R_table[i, j])){ # pick any; S and R should be filled in same slots
        # calculate effect of R logistic bacteria growth (abs) - only R growth if on antibiotics
        R_grow = r_growth * R_table[i-1, j] * (1 - ((R_table[i-1, j] + S_table[i-1, j])/total_capacity[i, j])) 
        # abx killing if abx.matrix is r abx (== 2) (abs)
        R_abx = -(abx.matrix[i-1, j] == 2) * abx.r * (R_table[i-1, j] + R_grow)
        # add effects to current table (abs first because log of a negative number is NaN)
        R_table[i, j] = R_table[i-1, j] + R_grow + R_abx
      } else {
        R_table[i, j] = R_table[i, j] #if filled (abs)
      }
      if (R_table[i, j] < 0) {R_table[i, j] = 0}
      
      if(is.na(S_table[i, j])){ 
        # calculate effect of S logistic bacteria growth (in absolute numbers)
        S_grow = s_growth * S_table[i-1, j] * (1 - ((S_table[i-1, j])/total_capacity[i, j])) # S may grow to max capacity
        # calculate killing effect of antibiotics R and antibiotics S (abs)
        S_abx = -(abx.matrix[i-1, j] >= 1) * abx.s * (S_table[i-1, j] + S_grow)
        # apply effects (abs first because log of a negative number is NaN)
        S_table[i, j] = S_table[i-1, j] + S_grow + S_abx
      } else { 
        S_table[i, j] = S_table[i, j] #if filled (abs)
      }
      if (S_table[i, j] < 0) {S_table[i, j] = 0}
      
      #transmission of R if roll pass prob check 
      roll = runif(1, 0, 1) # roll for transmission
      
      # transmission and trim range 
      ### transmission only happens if R and S have not exceeded total capacity and S is not fully occupying the capacity
      if(S_table[i, j] + R_table[i, j] > total_capacity[i, j]){ ## if existing S and R already exceed total capacity (abs)
        
        # if (S_table[i, j] > total_capacity[i, j]) { # S has a growing advantage if S alone exceeded capacity
        #   S_table[i, j] = total_capacity[i, j]
        #   R_table[i, j] = 0 
        # } else {
        #   S_table[i, j] = S_table[i, j]
        #   R_table[i, j] = total_capacity[i, j] - S_table[i, j] 
        # }
        
        S_table[i, j] =  S_table[i, j] / (S_table[i, j] + R_table[i, j]) * total_capacity[i, j]
        R_table[i, j] =  R_table[i, j] / (S_table[i, j] + R_table[i, j]) * total_capacity[i, j]

        #natural attrition of S and R take place according to their density if capacity exceeded
      
      } else { #transmission if S and R less than total_capacity
        
        R_trans = r_trans * R_table[i-1, j] * (roll < prop_r) # abs 
        
        if (S_table[i, j] + R_table[i, j] + R_trans > total_capacity[i, j]) {
          ### transmission if S and R less than total_capacity but with transmission exceeds total_capacity
          S_table[i, j] = S_table[i, j] / (S_table[i, j] + R_table[i, j] + R_trans) * total_capacity[i, j] #S and R+R_trans each given their proportion of the total capacity (abs)
          R_table[i, j] = (R_table[i, j] + R_trans) / (S_table[i, j] + R_table[i, j] + R_trans) * total_capacity[i, j]
          
        } else {
          ### transmission if sum of S and R and transmitted R are less than the total capacity
          #transmission happens with full amount of R_trans
          R_table[i, j] = R_table[i, j] + R_trans 
        }
      }
      
    }
  }
  
  total_capacity = log(total_capacity)
  r_thres_matrix = log(r_thres_matrix)
  S_table = log(S_table)
  R_table = log(R_table)
  
  return(list(S_table, R_table, total_capacity, r_thres_matrix))
}


