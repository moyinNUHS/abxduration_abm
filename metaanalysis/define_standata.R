# ========================================== #
# Define trial data for running JAGS models
# ========================================== #

dat = read.csv('lit_review/extracted_data.csv')
dat = dat[!is.na(dat$study_no),]

# dat = dat[-which(dat$abx_dur == 0),]
## function to write data into list for stan
write.jag.data <- function(fu_cutoff = 999, 
                           mean_prior,
                           mean_sd_prior,
                           sigma_prior1, 
                           sigma_prior2){
  
  dat = dat[which(dat$time_baseline_endoffu <= fu_cutoff),]
  
  dat$rat = dat$n_ind_outcome / dat$n_ind_contributedsamples
  plot(x = dat$abx_dur, y = dat$rat)
  
  ##############
  # Put into list 
  
  jags.data = list(
    pmid = dat$PMID,
    abx_dur = dat$abx_dur,
    time_baseline_endoffu = as.integer(dat$time_baseline_endoffu),
    n_ind_contributedsamples = round(dat$n_ind_contributedsamples), 
    n_ind_outcome = round(dat$n_ind_outcome), 
    study_id = dat$study_no, 
    study_N = max(dat$study_no), 
    setting = ifelse(dat$setting == 'Outpatient', 1, 0),
    mean_prior = mean_prior,
    mean_sd_prior = mean_sd_prior,
    sigma_prior1 = sigma_prior1, 
    sigma_prior2 = sigma_prior2
  )
  
  return(list(dat, jags.data))
  
}

