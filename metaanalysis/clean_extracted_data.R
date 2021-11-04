### clean extracted data from literature 

setwd('~/Documents/nBox/git_projects/abxduration_abm/metaanalysis/')

library(dplyr)
library(stringr)

dat = read.csv('lit_review/extracted_data.csv')

## remove leading and trailing spaces 
dat = as.data.frame(apply(dat, 2, str_trim))
dat = dat[!is.na(dat$study_no),] 

#############
# 5 trials #
unique(dat$study_no)
#############

############################################################################################
#### Get columns

dat_clean = dat[,c('PMID', #unique study identifier
                   'title',
                   'setting',
                   'abx_type',
                   'resistance_type',
                   'bacteria_reported',
                   'bacteria_grp',
                   'colonisation_site',
                   'abx_dur', #duration of antibiotics
                   'time_baseline_endoffu', # follow up time
                   'n_ind_contributedsamples',
                   'n_ind_outcome')]

write.csv(dat_clean, file = 'lit_review/clean_extracteddata.csv')
#############
# 7 trials #
#############

