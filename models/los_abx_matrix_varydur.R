# generate a table of number of days we want to observe (rows) -
# against number of beds in the ward (columns), filled in with patient id numbers

# the following function generates los and abx matrices when duration of antibiotics are varied, 
# number of courses of antibiotics between long and short roughly similar 

fill.abx.adm <- function(abx.type, meanDur, los){
  if (meanDur < los) {
    c(rep(abx.type, meanDur), rep(0, (los-meanDur)))
  } else {
    rep(abx.type, meanDur) # extend stay if duration exceeds los 
  }
}

los_abx_table_varydur <- function(n.bed, n.day, max.los, 
                                  p.infect, p.r.day1, p.r.after, p.infect.after, 
                                  meanDur, timestep, sd = 1){
  #define matrix structure 
  #rows - number of beds 
  #columns - days of observation 
  #assume all bed-days are occupied 
  n.bed = round(n.bed)
  n.day = round(n.day)
  max.beddays = n.bed * n.day
  
  ################
  # generate length of stay for each patient
  # (exponential distribution)
  ################
  all_los = ceiling(rexp(max.beddays, 1/max.los))
  max.patient = min(which(cumsum(all_los) > max.beddays)) + n.day
  all_los = all_los[1:max.patient]
  
  ################
  # generate antibiotic prescribing 
  # daily probability `p.infect` and`p.infect.after`
  ################
  patient.id = 1:max.patient
  adm.abx = lapply(all_los , function(x) {rep(0, x)})
  
  ## on admission 
  abx.adm.flip = rbinom(max.patient, 1, p.infect)
  abx.adm.id = patient.id[abx.adm.flip == 1] 
  abx.type = sample(c(2, 1), length(abx.adm.id), c(p.r.day1, (1 - p.r.day1)), replace = T) # p.r.day1 cannot be 0 
  abx.adm.id1 = abx.adm.id[abx.type == 1] # who got abx 1 
  abx.adm.id2 = abx.adm.id[abx.type == 2] # who got abx 2
  
  adm.abx[abx.adm.id2] = lapply(adm.abx[abx.adm.id2], function(x){
    fill.abx.adm(abx.type = 2, meanDur = meanDur, los = length(x))
  }) # fill abx 2
  
  adm.abx[abx.adm.id1] = lapply(adm.abx[abx.adm.id1], function(x){
    fill.abx.adm(abx.type = 1, meanDur = meanDur, los = length(x))
  }) # fill abx 1
  
  all_los = lengths(adm.abx)
  day1stay = c(1, cumsum(all_los) + 1)
  day1stay = day1stay[-length(day1stay)]
  
  los.expand = lapply(patient.id, function(x) {
    rep(x, all_los[x])}
  )
  
  ## during admission 
  abx.aft.flip = rbinom(sum(all_los), 1, p.infect.after)
  abx.aft.flip[day1stay] = 0 # change day of admission back to 0
  abx.aft.flip[abx.aft.flip > 0] = sample(c(2, 1), sum(abx.aft.flip), c(p.r.after, (1 - p.r.after)), replace = T)
  aft.abx = split(abx.aft.flip, unlist(los.expand)) 
  
  aft.abx = lapply(aft.abx, function(x)
    if (any(x > 0)){
      infect.days = which(x > 0)
      x.update = x 
      for (day1 in infect.days) {
        
        lastday = (day1 + meanDur - 1)
        x.update[day1 : lastday] = x[day1] # can be infected again while taking antibiotics
        
        #if (length(x) < lastday) break #if last day is beyond los, extend 
      }
      x.update
    } else {
      x
    }
  )
  
  ## fill adm.abx with 0 to the same los as aft.abx
  adm.abx = lapply(patient.id, function(i){
    length(adm.abx[[i]]) = length(aft.abx[[i]])
    adm.abx[[i]][is.na(adm.abx[[i]])] = 0
    return(adm.abx[[i]])
  })
  
  #update los
  all_los = lengths(adm.abx)
  
  ## make into matrix 
  patient.matrix = matrix(NA, nrow = n.day, ncol = n.bed)
  idx = 1
  for(j in 1:n.bed){

    sum_los = cumsum(all_los[idx:length(all_los)]) #cumulative stay of patients 
    
    #find last patient in the last row (last day of observation)
    end.pt.id = suppressWarnings(max(which(sum_los < n.day)) + idx)
    if (end.pt.id == -Inf) {end.pt.id = idx}
    if (end.pt.id > length(all_los)) {end.pt.id = length(all_los)}
    #fill matrix
    patient.matrix[, j] = rep(idx:end.pt.id, all_los[idx: end.pt.id])[1:n.day]
    
    idx = end.pt.id + 1
  }
  
  patient.id = 1:max(patient.matrix)
  los.expand = split(patient.matrix, unlist(patient.matrix))
  
  adm.abx = lapply(patient.id, function(i){
    length(adm.abx[[i]]) = length(los.expand[[i]])
    return(adm.abx[[i]])
  })
  aft.abx = lapply(patient.id, function(i){
    length(aft.abx[[i]]) = length(los.expand[[i]])
    return(aft.abx[[i]])
  })
  adm.abx = matrix(unlist(adm.abx), nrow = n.day, ncol = n.bed)
  aft.abx = matrix(unlist(aft.abx), nrow = n.day, ncol = n.bed)
  abx.matrix = pmax(adm.abx, aft.abx) # broad spectrum abx takes priority
  
  ################
  # find out if the patient was given antibiotics prior 
  # for 'treated individual outcome
  ################
  
  abx.b4 = split(abx.matrix, patient.matrix)
  abx.b4.fill = lapply(abx.b4, function(x){ 
    if (any(x >0)){
      x[min(which(x>0)):length(x)] = 1 
      x
    } else {
      x
    }
  })
  abxb4.matrix = matrix(unlist(abx.b4.fill), ncol = n.bed, nrow = n.day)
  
  return(list(patient.matrix, abx.matrix, abxb4.matrix))
}




