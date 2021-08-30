### clean extracted data from literature 

setwd('~/Documents/nBox/git_projects/indiv_abxduration_incmetaanalysis/indiv_abxduration/metaanalysis/')

library(dplyr)
library(stringr)

dat = read.csv('lit_review/extracted_data.csv')

## remove leading and trailing spaces 
dat = as.data.frame(apply(dat, 2, str_trim))
dat = dat[, c(1:which(colnames(dat) == 'n_ind_outcome'))]
dat = dat[!is.na(dat$study_no),] 

#############
# 26 trials #
#############

############################################################################################
###### group bacteria reported 
dat$bacteria_grp = NA
dat$bacteria_grp[which(dat$bacteria_reported == 'Acinectobacter baumannii' | 
                         dat$bacteria_reported == 'Pseudomonas' |
                         dat$bacteria_reported == 'Non fermenters')] = 'Non fermenters'
dat$bacteria_grp[which(dat$bacteria_reported == 'Citrobacter' | 
                         dat$bacteria_reported == 'E coli' |
                         dat$bacteria_reported == 'Edwardsiella' |
                         dat$bacteria_reported == 'Klebsiella pneumoniae' |
                         dat$bacteria_reported == 'Enterobacter' |
                         dat$bacteria_reported == 'Proteus' |
                         dat$bacteria_reported == 'Other Enterobacteriaceae' |
                         dat$bacteria_reported == 'Enterobacteriaceae' )] = 'Enterobacteriaceae'
dat$bacteria_grp[which(dat$bacteria_reported == 'Gram negatives' |
                         dat$bacteria_reported == 'H influenzae' |
                         dat$bacteria_reported == 'M catarhalis')] = 'Gram negatives'
dat$bacteria_grp[which(dat$bacteria_reported == 'Staph aureus'|
                         dat$bacteria_reported == 'Staph saprophyticus' |
                         dat$bacteria_reported == 'Pneumococcus' |
                         dat$bacteria_reported == 'Enterococcus')] = 'Gram positives'
dat$bacteria_grp[which(dat$bacteria_reported == 'Bacteria' | 
                         dat$bacteria_reported == 'Not reported')] = 'Unspecified'
dat$bacteria_grp[is.na(dat$bacteria_reported)] =  NA

############################################################################################
###### classify if the resistance detected could have been killed by the antibiotics used 
dat$abx_appropriate_for_resistancetype = NA

# if abx used covers the resistance type monitored = YES 
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Combination (cefotaxime, tobramycin, colistin, amphotericin B)' & 
                                               dat$resistance_type %in% c('Aminoglycoside', 'Fluoroquinolone', 'Cephalosporin', 
                                                                          'Resistance to >3 classes/antibiotics',  'HRMO',
                                                                          'Aminoglycosides, gentamicin, tobramycin, ciprofloxacin, ceftazidime',
                                                                          'Ceftazidime', 'Ciprofloxacin', 'Tobramycin', 'ESBL-producing'))] = 'YES'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Culture directed')] = 'YES'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Amox-clav' & dat$resistance_type == 'Penicillin')] = 'YES'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Azithromycin' & dat$resistance_type == 'Azithromycin')] = 'YES'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Combination (colistin, tobramycin, and nystatin)' & 
                                               dat$resistance_type %in% c('Resistance to >3 classes/antibiotics', 'ESBL-producing', 'HRMO', 'Gentamicin'))] = 'YES'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Tobramycin and ceftzidime' & 
                                               dat$resistance_type == 'Tobramycin')] = 'YES'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Nitrofurantoin, trimethoprim or cefalexin' & 
                                               dat$resistance_type %in% c('Cefalexin', 'Amox-clav', 'Nitrofurantoin', 'Trimethoprim',
                                                                          'Trimethoprim-sulfamethoxazole', 'Ciprofloxacin', 'Amoxicillin'))] = 'YES'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Combination (tobramycin and colistin sulphate)' & 
                                               dat$resistance_type %in% c('ESBL-producing', 'Tobramycin'))] = 'YES'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Combination (ciprofloxacin, gentamicin and polymyxin)' & 
                                               dat$resistance_type %in% c('Ciprofloxacin', 'Gentamicin'))] = 'YES'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Amoxicillin' & dat$resistance_type == 'Penicillin')] = 'YES'

# if abx used DOES NOT cover the resistance type monitored = NO 
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Combination (cefotaxime, tobramycin, colistin, amphotericin B)' & 
                                               dat$resistance_type %in% c('Imipenem', 'Polymyxin',
                                                                          'Vancomycin', 'Meticillin', 'Carbapenem', 'Colistin'))] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Culture directed' & 
                                               dat$resistance_type == 'Ceftazidime or imipenem-resistant')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Amoxicillin' & dat$resistance_type == 'Trimethoprim-sulfamethoxazole')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Azithromycin' & dat$resistance_type == 'Ampicillin')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Azithromycin' & dat$resistance_type == 'Meticillin')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Azithromycin' & dat$resistance_type == 'Penicillin')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Erythromycin' & 
                                               dat$resistance_type == 'Meticillin')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Ciprofloxacin' & dat$resistance_type == 'Fluoroquinolone')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Amox-clav' & dat$resistance_type == 'Beta-lactamase producing')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Combination (colistin, tobramycin, and nystatin)' & 
                                               dat$resistance_type %in% c('Carbapenem', 
                                                                          'Colistin', 'Vancomycin', 'Meticillin'))] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Trimethoprim-sulfamethoxazole' & dat$resistance_type == 'Trimethoprim-sulfamethoxazole')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Moxfloxacin' & 
                                               dat$resistance_type == 'Fluoroquinolone')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Amoxicillin' & dat$resistance_type == 'Beta-lactamase producing')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Trimethoprim-sulfamethoxazole' & dat$resistance_type == 'Ceftibuten')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Ceftibuten' & dat$resistance_type == 'Trimethoprim-sulfamethoxazole')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Ceftibuten' & dat$resistance_type == 'Ceftibuten')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Fosfomycin trometamol' & dat$resistance_type == 'Fosfomycin trometamol')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Ciprofloxacin' & dat$resistance_type == 'Ciprofloxacin')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Nitrofurantoin, trimethoprim or cefalexin' & 
                                               dat$resistance_type == 'Mecillinam')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Amoxicillin' & dat$resistance_type == 'ESBL-producing')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Flucloxacillin' & 
                                               dat$resistance_type == 'ESBL-producing')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Flucloxacillin' & 
                                               dat$resistance_type == 'Meticillin')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Ceftiaxone' & dat$resistance_type == 'ESBL-producing')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Ceftiaxone' & dat$resistance_type == 'Meticillin')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Fosfomycin trometamol or trimethoprim-sulfamethoxazole or amox-clav or cefixim or furadantin' & 
                                               dat$resistance_type == 'ESBL-producing')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Fosfomycin trometamol or trimethoprim-sulfamethoxazole or amox-clav or cefixim or furadantin' & 
                                               dat$resistance_type == 'Meticillin')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Meropenem' & 
                                               dat$resistance_type == 'Carbapenem')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Ampicillin+gentamicin or cefotaxime+gentamicin' & dat$resistance_type == 'Carbapenem')] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Combination (tobramycin and colistin sulphate)' & 
                                               dat$resistance_type %in% c('Carbapenem', 'Colistin'))] = 'NO'
dat$abx_appropriate_for_resistancetype[which(dat$abx_type == 'Combination (ciprofloxacin, gentamicin and polymyxin)' & 
                                               dat$resistance_type %in% c('Oxacillin', 'Polymyxin'))] = 'NO'


#dat[which(is.na(dat$abx_appropriate_for_resistancetype) & !is.na(dat$abx_type)),]


############################################################################################
#### classify if the antibiotic resistance type has effective treatment or not 
dat$abx_appropriate_exist = 'YES' 
dat$abx_appropriate_exist[which(dat$resistance_type %in% c('Ceftazidime or imipenem-resistant', 
                                                              'Imipenem', 'Carbapenem',
                                                              'Colistin', 'Polymyxin'))] = 'NO'

#############
length(unique(dat$PMID))
# 24 trials #
#############

############################################################################################
#### decide which type of study 
# 1. those with clear follow up period (days between antibiotic intake and surveillance sample) - binomial distribution
# 2. those with unclear follow up period (point prevalence) - poisson distribution 

dat$outcome_distribution = 'binomial'
dat$outcome_distribution[which(dat$time_baseline_endoffu == 'Point prevalence surveys')] = 'poisson'

############################################################################################
#### decide which type of study 

############################################################################################
#### Get columns

dat_clean = dat[,c('PMID', #unique study identifier
                   'title',
                   'setting',
                   'abx_appropriate_for_resistancetype', #separate into scenarios where appropriate antibiotics were used or not
                   'abx_appropriate_exist',
                   'abx_type',
                   'resistance_type',
                   'bacteria_reported',
                   'bacteria_grp',
                   'colonisation_site',
                   'abx_dur', #duration of antibiotics
                   'time_baseline_endoffu', # follow up time
                   'n_ind_contributedsamples',
                   'n_ind_outcome',
                   'outcome_distribution')] # distribution to use to model the outcome

write.csv(dat_clean, file = 'lit_review/clean_extracteddata.csv')
#############
# 24 trials #
#############

