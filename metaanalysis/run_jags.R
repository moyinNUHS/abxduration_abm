rm(list = ls())

library(R2jags)
library(matrixStats)
library(loo)
library(ggmcmc)

source('define_standata.R')
source('jags.R')

##########################
### No appropriate antibiotics 
##########################

# MAIN ANALYSIS WITH 30 days followup cut off, all settings #

#get data 
out = select.data(
  abx_appropriate_for_resistancetype = 'NO',
  bacteria_grp = c('Enterobacteriaceae', 'Non fermenters', 'Gram negatives')
)

dftocheck = out[[1]]
lapply(unique(dftocheck$PMID), function(x){
  unique(dftocheck $time_baseline_endoffu[which(dftocheck$PMID == x)])
})

jags.data = out[[2]]
print(paste('There are', max(lengths(jags.data)), 'datapoints from', jags.data$study_N, 'trials.'))

# run model 
fit <- jags(data = jags.data, inits = init, parameters.to.save = params, model.file = jags.mod,
            n.chains = 4, n.iter = 500000, n.burnin = 100000, n.thin = 100, DIC = T)

jags.out = list(jags.data, fit)

save(jags.out, file = paste0('runs/jagsOUT_noappro_30d_all_', Sys.Date(), '.Rdata'))


# SENSITIVITY ANALYSIS WITH 60 days followup cut off, all settings #
out = select.data(
  abx_appropriate_for_resistancetype = 'NO',
  fu_cutoff = 60,
  bacteria_grp = c('Enterobacteriaceae', 'Non fermenters', 'Gram negatives')
)

dftocheck = out[[1]]
jags.data = out[[2]]
print(paste('There are', max(lengths(jags.data)), 'datapoints from', jags.data$study_N, 'trials.'))
fit <- jags(data = jags.data, inits = init, parameters.to.save = params, model.file = jags.mod,
            n.chains = 4, n.iter = 500000, n.burnin = 100000, n.thin = 100, DIC = T)

jags.out = list(jags.data, fit)

save(jags.out, file = paste0('runs/jagsOUT_noappro_60d_all_', Sys.Date(), '.Rdata'))

# SENSITIVITY ANALYSIS WITH 30 days followup cut off, hosp settings #
out = select.data(
  abx_appropriate_for_resistancetype = 'NO',
  setting = c("ICU", "Neonatal unit", "Hospital"),
  bacteria_grp = c('Enterobacteriaceae', 'Non fermenters', 'Gram negatives')
)

dftocheck = out[[1]]
jags.data = out[[2]]
print(paste('There are', max(lengths(jags.data)), 'datapoints from', jags.data$study_N, 'trials.'))
fit <- jags(data = jags.data, inits = init, parameters.to.save = params, model.file = jags.mod,
            n.chains = 4, n.iter = 500000, n.burnin = 100000, n.thin = 100, DIC = T)

jags.out = list(jags.data, fit)

save(jags.out, file = paste0('runs/jagsOUT_noappro_30d_hosp_', Sys.Date(), '.Rdata'))


###############################################
### With appropriate antibiotics 
###############################################

out = select.data(
  #abx_appropriate_exist = 'NO', 
  abx_appropriate_for_resistancetype = 'YES',
  #colonisation_site = c('Gut', 'Gut and respiratory tract',  'Gut or respiratory tract', 'Urine'), 
  bacteria_grp = c('Enterobacteriaceae', 'Non fermenters', 'Gram negatives')
)

dftocheck = out[[1]]
lapply(unique(dftocheck$PMID), function(x){
  unique(dftocheck$time_baseline_endoffu[which(dftocheck$PMID == x)])
})

jags.data = out[[2]]
print(paste('There are', max(lengths(jags.data)), 'datapoints from', jags.data$study_N, 'trials.'))

fit <- jags(data = jags.data, inits = init, parameters.to.save = params, model.file = jags.mod,
            n.chains = 4, n.iter = 500000, n.burnin = 100000, n.thin = 100, DIC = T)

jags.out = list(jags.data, fit)

save(jags.out, file = paste0('runs/jagsOUT_appro_30d_all_', Sys.Date(), '.Rdata'))

# SENSITIVITY ANALYSIS WITH 60 days followup cut off, all settings #
out = select.data(
  abx_appropriate_for_resistancetype = 'YES',
  fu_cutoff = 60,
  bacteria_grp = c('Enterobacteriaceae', 'Non fermenters', 'Gram negatives')
)

dftocheck = out[[1]]
jags.data = out[[2]]
print(paste('There are', max(lengths(jags.data)), 'datapoints from', jags.data$study_N, 'trials.'))
fit <- jags(data = jags.data, inits = init, parameters.to.save = params, model.file = jags.mod,
            n.chains = 4, n.iter = 500000, n.burnin = 100000, n.thin = 100, DIC = T)

jags.out = list(jags.data, fit)

save(jags.out, file = paste0('runs/jagsOUT_appro_60d_all_', Sys.Date(), '.Rdata'))

# SENSITIVITY ANALYSIS WITH 30 days followup cut off, hosp settings #
out = select.data(
  abx_appropriate_for_resistancetype = 'YES',
  setting = c("ICU", "Neonatal unit", "Hospital"),
  bacteria_grp = c('Enterobacteriaceae', 'Non fermenters', 'Gram negatives')
)

dftocheck = out[[1]]
jags.data = out[[2]]
print(paste('There are', max(lengths(jags.data)), 'datapoints from', jags.data$study_N, 'trials.'))
fit <- jags(data = jags.data, inits = init, parameters.to.save = params, model.file = jags.mod,
            n.chains = 4, n.iter = 500000, n.burnin = 100000, n.thin = 100, DIC = T)

jags.out = list(jags.data, fit)

save(jags.out, file = paste0('runs/jagsOUT_appro_30d_hosp_', Sys.Date(), '.Rdata'))

###############################################
### check models 
###############################################

###1. View trace plots
traceplot(fit, varname = c('a'))
traceplot(fit_PT_noint_ward_main[[1]], varname = c('b'))
traceplot(fit_PT_noint_ward_main[[1]], varname = c('c'))
traceplot(fit_PT_noint_ward_main[[1]], varname = c('d'))

traceplot(fit_PT_noint_ward_main[[1]], varname = c('a'))
traceplot(fit_PT_noint_ward_main[[1]], varname = c('b'))
traceplot(fit_PT_noint_ward_main[[1]], varname = c('c'))
traceplot(fit_PT_noint_ward_main[[1]], varname = c('d'))

###2. Check R hat and effective chains
summary(fit$BUGSoutput$summary)
summary(fit_STAFF_noint_ward_main[[1]]$BUGSoutput$summary)

