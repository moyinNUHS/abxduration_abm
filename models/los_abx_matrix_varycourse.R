# generate a table of number of days we want to observe (rows) -
# against number of beds in the ward (columns), filled in with patient id numbers
# accommodating for increase in length of stay with hospital acquired infections

# the following function generates los and abx matrices when the number of courses of antibiotics varied,
# duration of antibiotics are the same 

# differences with los_abx_matrix_varydur: 
# 1. define abx.type (instead of using p.r.after) - such that interrupted courses are for 1 type of antibiotics only 
# 2. generates 2 abx matrixes - one with continuous course and one with interrupted courses - since los is the same
#    in both arms

load("models/norm_ecdf.Rdata")

###################
# define los and antibiotic duration  #####################
los_abx_table_course <- function(n.bed, n.day, max.los, abx.type, 
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

