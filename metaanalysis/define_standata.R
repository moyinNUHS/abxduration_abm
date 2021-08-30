# ========================================== #
# Define trial data for running JAGS models
# ========================================== #

dat = read.csv('metaanalysis/lit_review/clean_extracteddata.csv')

## function to write data into list for stan
select.data <- function(abx_appropriate_for_resistancetype = unique(dat$abx_appropriate_for_resistancetype), 
                        abx_appropriate_exist = unique(dat$abx_appropriate_exist), 
                        colonisation_site = unique(dat$colonisation_site), 
                        setting = unique(dat$setting),
                        fu_cutoff = 30, #main analysis use 30 days as cut off for follow up time 
                        bacteria_grp = unique(dat$bacteria_grp)){
  
  ##############
  # Select data from main dataset 
  #### select rows for columns with no NA 
  dat.to.analyse = dat[which(dat$colonisation_site %in% colonisation_site & 
                               dat$bacteria_grp %in% bacteria_grp & 
                               dat$setting %in% setting &
                               dat$abx_appropriate_exist %in% abx_appropriate_exist),]
  
  #### select rows based on abx_appropriate_for_resistancetype which includes NA values for dur = 0 
  dat.to.analyse = dat.to.analyse[- which(dat.to.analyse$abx_appropriate_for_resistancetype != abx_appropriate_for_resistancetype),]
  # hence only left with rows with NA and value for abx_appropriate_for_resistancetype
  
  #### select rows depending on follow up time 
  fu.time = NA
  fu.time[which(dat.to.analyse$time_baseline_endoffu != 'Point prevalence surveys')] = as.numeric(dat.to.analyse$time_baseline_endoffu[which(dat.to.analyse$time_baseline_endoffu != 'Point prevalence surveys')])
  fu.time[which(dat.to.analyse$time_baseline_endoffu == 'Point prevalence surveys')] = 1
  dat.to.analyse = dat.to.analyse[which(fu.time <= fu_cutoff),]
  
  #### remove trials with just one duration 
  toremove = unlist(lapply(unique(dat.to.analyse$PMID), function(x){
    abx.dur = unique(dat.to.analyse$abx_dur[which(dat.to.analyse$PMID == x)])
    if(length(abx.dur) == 1){x}
  }))
  
  if (length(toremove) > 0) {dat.to.analyse = dat.to.analyse[-which(dat.to.analyse$PMID %in% toremove),]}
  
  # then select the matching rows with !is.na(dat.to.analyse$abx_appropriate_for_resistancetype) that has abx.dur = 0
  rowstokeep = list()
  dat.to.analyse$PMID = as.character(dat.to.analyse$PMID)
  for (s in unique(dat.to.analyse$PMID)){
    
    res = unique(dat.to.analyse[which(dat.to.analyse$PMID == s & !is.na(dat.to.analyse$abx_appropriate_for_resistancetype)), 'resistance_type'])
    bact = unique(dat.to.analyse[which(dat.to.analyse$PMID == s & !is.na(dat.to.analyse$abx_appropriate_for_resistancetype)),  'bacteria_grp'])
    rowstokeep[[s]] = which(dat.to.analyse$PMID == s & 
                              dat.to.analyse$resistance_type %in% res & 
                              dat.to.analyse$bacteria_grp %in% bact)
    
  }
  dat.to.analyse = dat.to.analyse[unlist(rowstokeep),]
  
  dat.to.analyse$study_id = cumsum(c(0, diff(as.numeric(dat.to.analyse$PMID))) != 0) + 1
  
  dat.to.analyse$rat = dat.to.analyse$n_ind_outcome/ dat.to.analyse$n_ind_contributedsamples
  plot(x = dat.to.analyse$abx_dur, y = dat.to.analyse$rat)
  
  ##############
  # Put into list 
  
  jags.data = list(
    pmid = dat.to.analyse$PMID,
    abx_dur = dat.to.analyse$abx_dur,
    time_baseline_endoffu = ifelse(dat.to.analyse$time_baseline_endoffu == "Point prevalence surveys", -99, as.integer(dat.to.analyse$time_baseline_endoffu)),
    n_ind_contributedsamples = round(dat.to.analyse$n_ind_contributedsamples), 
    n_ind_outcome = round(dat.to.analyse$n_ind_outcome), 
    binom_dist = which(dat.to.analyse$outcome_distribution == 'binomial'),
    pois_dist = which(dat.to.analyse$outcome_distribution == 'poisson'),
    study_id = dat.to.analyse$study_id, 
    study_N = max(dat.to.analyse$study_id)
  )
  
  
  return(list(dat.to.analyse, jags.data))
  
}

