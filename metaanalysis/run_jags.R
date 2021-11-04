rm(list = ls())

library(R2jags)
library(matrixStats)
library(loo)
library(ggmcmc)

setwd('~/Documents/nBox/git_projects/abxduration_abm/metaanalysis/')
source('define_standata.R')
source('jags.R')

#get data 
out = write.jag.data (mean_prior = 0, 
                      mean_sd_prior = 0.1, # sqrt(1/0.1) = 3 
                      sigma_prior1 = 1,  
                      sigma_prior2 = 0.1
)

dftocheck = out[[1]]
jags.data = out[[2]]
print(paste('There are', max(lengths(jags.data)), 'datapoints from', jags.data$study_N, 'trials.'))

# run model 

##### WITH RANDOM EFFECT AND D
fit <- jags(data = jags.data, inits = init.rand, parameters.to.save = params.rand, model.file = jags.rand,
            n.chains = 4, n.iter = 800000, n.burnin = 100000, n.thin = 100, DIC = T)
jags.out = list(jags.data, fit)
save(jags.out, file = paste0('runs/rand_0313_', Sys.Date(), '.Rdata'))

##### WITH RANDOM EFFECT AND NO D (setting)
fit <- jags(data = jags.data, inits = init.randnod, parameters.to.save = params.randnod, model.file = jags.randnod,
            n.chains = 4, n.iter = 800000, n.burnin = 100000, n.thin = 100, DIC = T)
jags.out = list(jags.data, fit)
save(jags.out, file = paste0('runs/randnod_0313_', Sys.Date(), '.Rdata'))

##### NO RANDOM EFFECTS
fit <- jags(data = jags.data, inits = init.norand, parameters.to.save = params.norand, model.file = jags.norand,
            n.chains = 4, n.iter = 800000, n.burnin = 100000, n.thin = 100, DIC = T)
jags.out = list(jags.data, fit)
save(jags.out, file = paste0('runs/norand_0313_', Sys.Date(), '.Rdata'))

##### CHANGE PRIOR
out = write.jag.data(
  mean_prior = 0, 
  mean_sd_prior = 0.02,   # sqrt(1/0.1) = 7
  sigma_prior1 = 1, 
  sigma_prior2 = 0.02
)
jags.data = out[[2]]
fit <- jags(data = jags.data, inits = init.randnod, parameters.to.save = params.randnod, model.file =  jags.randnod,
            n.chains = 4, n.iter = 800000, n.burnin = 100000, n.thin = 100, DIC = T)
jags.out = list(jags.data, fit)
save(jags.out, file = paste0('runs/randnod_0717_', Sys.Date(), '.Rdata'))

###### Sensitivity analysis with 30 days cut off 
#get data 
out = write.jag.data(
  fu_cutoff = 30, 
  mean_prior = 0, 
  mean_sd_prior = 0.1, # sqrt(1/0.1) = 3 
  sigma_prior1 = 1,  
  sigma_prior2 = 0.1 
)
jags.data = out[[2]]
fit <- jags(data = jags.data, inits = init.randnod, parameters.to.save = params.randnod, model.file = jags.randnod,
            n.chains = 4, n.iter = 800000, n.burnin = 100000, n.thin = 100, DIC = T)
jags.out = list(jags.data, fit)
save(jags.out, file = paste0('runs/randnodsens_0313_', Sys.Date(), '.Rdata'))


###############################################
### check models 
###############################################
load('runs/rand_0313_2021-11-02.Rdata')
fit = jags.out[[2]]

###1. View trace plots
traceplot(fit, varname = c('a'))
traceplot(fit, varname = c('b'))
traceplot(fit, varname = c('c'))
traceplot(fit, varname = c('d'))

par(ask=F) 

###2. Check R hat and effective chains
summary(fit$BUGSoutput$summary)

###3. Check WAIC 
# prepare data 
check_waic = list()
files = list.files('runs')[grep('2021-11-02', list.files('runs'))]
for (file in files){
  
  load(paste0('runs/', file))
  
  all_matrix = jags.out[[2]]$BUGSoutput$sims.matrix
  loglik_matrix_list = all_matrix[,grep("loglike", colnames(all_matrix))]
  waic = t(as.data.frame(waic(loglik_matrix_list, pointwise = TRUE)$estimates))['Estimate','waic']
  dic = jags.out[[2]]$BUGSoutput$DIC
  
  check_waic[[file]] = c(model = file, waic = round(waic, 2), dic = round(dic, 2))
}

compare = do.call('rbind', check_waic)
compare[1:3,]
