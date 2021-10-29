rm(list = ls())

library(R2jags)
library(matrixStats)
library(loo)
library(ggmcmc)

setwd('~/Documents/nBox/git_projects/abxduration_abm/metaanalysis/')
source('define_standata.R')
source('jags.R')

#get data 
out = select.data(
  abx_appropriate_for_resistancetype = 'NO',
  fu_cutoff = 60,
  bacteria_grp = c('Enterobacteriaceae', 'Gram negatives'),
  mean_prior = 0, 
  mean_sd_prior = 0.1, # sqrt(1/0.1) = 3 
  sigma_prior1 = 1,  
  sigma_prior2 = 0.1 
)

dftocheck = out[[1]]
jags.data = out[[2]]
print(paste('There are', max(lengths(jags.data)), 'datapoints from', jags.data$study_N, 'trials.'))

# run model 

##### MAIN MODEL 
fit <- jags(data = jags.data, inits = init.rand, parameters.to.save = params.rand, model.file = jags.rand,
            n.chains = 4, n.iter = 500000, n.burnin = 100000, n.thin = 100, DIC = T)
jags.out = list(jags.data, fit)
save(jags.out, file = paste0('runs/rand_0313_', Sys.Date(), '.Rdata'))

##### NO RANDOM EFFECTS
fit <- jags(data = jags.data, inits = init.norand, parameters.to.save = params.norand, model.file = jags.norand,
            n.chains = 4, n.iter = 500000, n.burnin = 100000, n.thin = 100, DIC = T)
jags.out = list(jags.data, fit)
save(jags.out, file = paste0('runs/norand_0313_', Sys.Date(), '.Rdata'))

##### NO D (setting)
fit <- jags(data = jags.data, inits = init.norandnod, parameters.to.save = params.norandnod, model.file = jags.norandnod,
            n.chains = 4, n.iter = 500000, n.burnin = 100000, n.thin = 100, DIC = T)
jags.out = list(jags.data, fit)
save(jags.out, file = paste0('runs/norand_0313nod_', Sys.Date(), '.Rdata'))

##### CHANGE PRIOR
out = select.data(
  abx_appropriate_for_resistancetype = 'NO',
  fu_cutoff = 60,
  bacteria_grp = c('Enterobacteriaceae', 'Gram negatives'),
  mean_prior = 0, 
  mean_sd_prior = 0.02,   # sqrt(1/0.1) = 7
  sigma_prior1 = 1, 
  sigma_prior2 = 0.02
)
jags.data = out[[2]]
fit <- jags(data = jags.data, inits = init.norandnod, parameters.to.save = params.norandnod, model.file =  jags.norandnod,
            n.chains = 4, n.iter = 500000, n.burnin = 100000, n.thin = 100, DIC = T)
jags.out = list(jags.data, fit)
save(jags.out, file = paste0('runs/norand_0717nod_', Sys.Date(), '.Rdata'))

###### Sensitivity analysis with 30 days cut off 
#get data 
out = select.data(
  abx_appropriate_for_resistancetype = 'NO',
  fu_cutoff = 30,
  bacteria_grp = c('Enterobacteriaceae', 'Gram negatives'),
  mean_prior = 0, 
  mean_sd_prior = 0.1, # sqrt(1/0.1) = 3 
  sigma_prior1 = 1,  
  sigma_prior2 = 0.1 
)
jags.data = out[[2]]
fit <- jags(data = jags.data, inits = init.norandnod, parameters.to.save = params.norandnod, model.file = jags.norandnod,
            n.chains = 4, n.iter = 500000, n.burnin = 100000, n.thin = 100, DIC = T)
jags.out = list(jags.data, fit)
save(jags.out, file = paste0('runs/norandsens_0313nod_', Sys.Date(), '.Rdata'))


###############################################
### check models 
###############################################
load('runs/norand_0313_2021-10-27.Rdata')
fit = jags.out[[2]]

###1. View trace plots
traceplot(fit, varname = c('a'))
traceplot(fit, varname = c('b'))
traceplot(fit, varname = c('c'))

par(ask=F) 

###2. Check R hat and effective chains
summary(fit$BUGSoutput$summary)

###3. Check WAIC 
# prepare data 
check_waic = list()
files = list.files('runs')[grep('rand_0313', list.files('runs'))]
for (file in files){
  
  load(paste0('runs/', file))
  
  all_matrix = jags.out[[2]]$BUGSoutput$sims.matrix
  loglik_matrix_list = all_matrix[,grep("loglike", colnames(all_matrix))]
  
  check_waic[[file]] = waic(loglik_matrix_list, pointwise = TRUE)
}

out.waic = cbind(check_waic[[1]]$estimates, check_waic[[2]]$estimates, check_waic[[3]]$estimates)
colnames(out.waic) = c('No random effect est', 'No random effect SE', 
                       'No random effect noD est', 'No random effect noD SE', 
                       'With random effect estnod', 'With random effect SEnod')
out.waic
