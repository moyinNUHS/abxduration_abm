# generate a table of number of days we want to observe (rows) -
# against number of beds in the ward (columns), filled in with patient id numbers
# accommodating for increase in length of stay with hospital acquired infections

# the following function generates los and abx matrices when duration of antibiotics are varied, 
# number of courses of antibiotics between long and short roughly similar 

# differences with los_abx_matrix_varycourse: 
# 1. define p.r.after (instead of using abx.type ) - such that antibiotics prescribed after admission may vary 
# 2. generates 1 abx matrix only -since los will differ in both arms

load("~/Documents/nBox/git_projects/indiv_abxduration_incmetaanalysis/indiv_abxduration/models/norm_ecdf.Rdata")

los_abx_table_varydur <- function(n.bed, n.day, max.los, 
                                  p.infect, p.r.day1, p.r.after, cum.r.1, 
                                  meanDur, timestep, sd = 1){
  
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
  no.abx.r.day1 = round(n.patient * p.infect * p.r.day1) #number of broad spectrum spectrum antibiotics prescribed for community acquired infection 
  no.abx.s.day1 = round(n.patient * p.infect) - no.abx.r.day1 #number of narrow spectrum antibiotics prescribed for community acquired infection 
  #identify patient ids who were prescribed broad vs narrow antibiotics
  id.abx.r.day1 = sample(patient.id, no.abx.r.day1) #id of patients prescribed broad spectrum antibiotics on day 1
  id.abx.s.day1 = sample(patient.id[-id.abx.r.day1], no.abx.s.day1) #id of patients prescribed narrow spectrum antibiotics on day 1
  
  # then decide if going to receive antibiotics after admission - risk of acquiring HAI increases with length of stay 
  # (cumulative normal distribution)
  # if prescribed antibiotics after admission then prolong admission to complete course 
  prob.after = norm_ecdf((1/cum.r.1) * (1 : max(all_los))) # get the cumulative probability for each day
  # norm_ecdf - imported normalized function
  abx.dayafter = lapply(patient.id, function(id){
    # if the patient will receive antibiotics after admission - 1 if will be prescribed antibiotics on the day during los 
    abx.dayafter.binary = rbinom(all_los[id], 1, prob = prob.after[1:all_los[id]])
    # start dates if given antibiotics
    if (sum(abx.dayafter.binary) > 0) {which(abx.dayafter.binary == 1)} else { 0 }
  })
  id.abx.dayafter = patient.id[which(lengths(abx.dayafter) > 1)]
  
  #allocate antibiotics for the patients upon and during admission
  all_abx = lapply(patient.id, function(id){ # for each patient 
    
    # create empty vector with 0 with length same as length of stay
    fill = rep(0, all_los[id])
    
    # fill antibiotics on admission
    if (length(id.abx.s.day1) > 0 & id %in% id.abx.s.day1) { # if the patient is given narrow spectrum antibiotics on admission
      # antibiotic duration randomly drawn from a truncated normal distribution
      dur = round(rtnorm(1, mean = meanDur, lower = 1, sd = sd))
      fill[1: dur] = 1
    } else if (length(id.abx.r.day1) > 0 & id %in% id.abx.r.day1) {
      # antibiotic duration randomly drawn from a truncated normal distribution
      dur = round(rtnorm(1, mean = meanDur, lower = 1, sd = sd))
      fill[1: dur] = 2
    }
    
    # fill antibiotics after admission
    if (length(id.abx.dayafter) > 0 & id %in% id.abx.dayafter) { # if the patient is given antibiotics after admission
      
      # how many courses and duration for each course 
      courses.start = abx.dayafter[[id]][which(abx.dayafter[[id]] > 0)]
      no.courses = length(courses.start)
      dur = round(rtnorm(no.courses, mean = meanDur, lower = 1, sd = sd))
      
      # which antibiotics for which course 
      no.r.courses = round(no.courses * p.r.after)
      r.courses.start = sample(courses.start, no.r.courses)
      
      for (start in courses.start){
        narrow.or.broad = ifelse(start %in% r.courses.start, 2, 1)
        end = start + dur[which(start == courses.start)] - 1
        fill[start : end] = narrow.or.broad 
      }
      
    }
    
    return(fill)
  })
  
  ####################
  #PUT ABX COURSES AND LOS INTO PATIENT and ABX MATRIX
  #make up a matrix of number of days we want to observe (rows) -
  #against number of beds in the ward (columns)
  all_los = lengths(all_abx)#update length of stay incorporating abx.r
  patient.matrix = abx.matrix = matrix(NA, nrow = n.day, ncol = n.bed)
  idx = 1
  for(j in 1:n.bed){
    
    sum_los = cumsum(all_los[idx:length(all_los)]) #cumulative stay of patients 
    
    #find last patient in the last row (last day of observation)
    end.pt.id = suppressWarnings(max(which(sum_los < n.day)) + idx)
    if (end.pt.id == -Inf) {end.pt.id = idx}
    if (end.pt.id > length(all_los)) {end.pt.id = length(all_los)}
    #fill matrix
    patient.matrix[, j] = rep(idx:end.pt.id, all_los[idx: end.pt.id])[1:n.day]
    abx.matrix[,j] = unlist(all_abx[idx:end.pt.id])[1:n.day]
    
    idx = end.pt.id + 1
  }
  
  patient.matrix = patient.matrix[rep(1:n.day, each = timestep), ]
  abx.matrix = abx.matrix[rep(1:n.day, each = timestep), ]
  
  return(list(patient.matrix, abx.matrix))
}




