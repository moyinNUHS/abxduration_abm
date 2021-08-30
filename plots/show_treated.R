# plot R with days amongst the treated 

library(ggplot2)

## load data
dat.names.list = list('runs/simple3state_500_treated_zero_2021-07-01.Rdata',
                      'runs/simple3state_500_treated_notzero_2021-07-02.Rdata',
                      'runs/cocarriage5state_500_treated_zero_2021-07-02.Rdata',
                      'runs/cocarriage5state_500_treated_notzero_2021-07-01.Rdata',
                      'runs/populationgrowth_500_treated_zero_2021-06-30.Rdata',
                      'runs/populationgrowth_500_treated_notzero_2021-06-30.Rdata')
lapply(dat.names.list[grep('_zero', dat.names.list)], load, .GlobalEnv)
LHS.simple.zero = LHS.simple
LHS.cocarriage.zero = LHS.cocarriage
LHS.populationgrowth.zero = LHS.populationgrowth
lapply(dat.names.list[grep('notzero', dat.names.list)], load, .GlobalEnv)
LHS.simple.notzero = LHS.simple
LHS.cocarriage.notzero = LHS.cocarriage
LHS.populationgrowth.notzero = LHS.populationgrowth

d.list = list(LHS.simple.zero, LHS.cocarriage.zero, LHS.populationgrowth.zero,
              LHS.simple.notzero, LHS.cocarriage.notzero, LHS.populationgrowth.notzero)
rm(list = ls()[-which(ls() == 'd.list')])

## clean data 

coltoretrieve.abxr.pres = c('p.r.day1','p.r.after')
coltoretrieve.transmission = c('pi_ssr')
coltoretrieve.baselineR = c('prop_R')
coltoretrieve.bactgrowth = c('fitness.r', 's_growth', 'r_growth', 'repop.s')

intersect_all <- function(list){
  Reduce(intersect, list)
}

classify <- function(dat, colnames, type) {
  
  coltoretrieve = colnames(dat)[colnames(dat) %in% colnames]
  cols = dat[,c(coltoretrieve)]
  if (length(coltoretrieve) > 1) {
    
    cutoff = lapply(cols, quantile, c(0.1, 0.9))
    less = list()
    more = list()
    for (i in 1:length(coltoretrieve)){
      less[[i]] = which(cols[,i] < cutoff[[i]][1])
      more[[i]] = which(cols[,i] > cutoff[[i]][2])
    }
    
    out = rep(NA, nrow(dat))
    out[intersect_all(less)] = paste0('low', type)
    out[intersect_all(more)] = paste0('high', type)
    
  } else {
    
    cutoff = quantile(cols, c(0.2, 0.8))
    less = which(cols < cutoff[1])
    more = which(cols > cutoff[2])
    out = rep(NA, nrow(dat))
    out[less] = paste0('low', type)
    out[more] = paste0('high', type)
  }
  

  
  return(out)
  
}

clean.d.list = lapply(d.list, function(x){
  
  abxr.pres = classify(dat = x$data, colnames = coltoretrieve.abxr.pres, type = 'abxr.pres')
  transmission = classify(dat = x$data, colnames = coltoretrieve.transmission, type = 'transmission')
  baselineR = classify(dat = x$data, colnames = coltoretrieve.baselineR, type = 'baselineR')
  bactgrowth = classify(dat = x$data, colnames = coltoretrieve.bactgrowth, type = 'bactgrowth')
  
  y = x$res[,1:20,1]
  y = cbind(rep(0, nrow(y)), y)

  data.frame(abxr.pres = rep(abxr.pres, each = 21), 
             transmission = rep(transmission, each = 21), 
             baselineR = rep(baselineR, each = 21), 
             bactgrowth = rep(bactgrowth, each = 21), 
             day = rep(0:20, each = nrow(y)),
             model = as.character(x$call[[2]]),
             scenario = ifelse(max(x$data$abx.r < 0.1), 
                               'Resistance phenotype treated with ineffective antibiotic', 
                               'Resistance phenotype treated with effective antibiotic'),
             risk = as.vector(y))
  
})
plot.d = do.call('rbind', clean.d.list)

plot.d$scenario = as.factor(plot.d$scenario)
plot.d$scenario = factor(plot.d$scenario, 
                         levels = c('Resistance phenotype treated with ineffective antibiotic', 
                                    'Resistance phenotype treated with effective antibiotic'))

p1 = ggplot(data = plot.d, aes(x = day, y = risk)) + 
  geom_point(data = plot.d, aes(color = model), alpha = 0.3) + 
  geom_smooth(color = 'grey20', se = F) +  
  labs(y = 'Probability of colonisation by resistant Gram negative bacteria\ncompared to baseline prevalence of resistance carriers\namongst the newly admitted patients', 
       x = 'Antibiotic duration (days)') +
  scale_color_manual(values = c('#FAC748', '#90A9B7', '#E55934'), name = 'Models', labels = c('Simple 3-state', 'Co-carriage 5 state', 'Population growth')) + 
  facet_wrap(~ scenario) + 
  theme_minimal() + 
  theme(legend.position = 'bottom', 
        text = element_text(size = 15))
p1 

ggplot(data = plot.d) + 
  geom_point(data = subset(plot.d, !is.na(bactgrowth)), aes(x = day, y = risk, color = bactgrowth), alpha = 0.3) + 
  geom_smooth(data = subset(plot.d, !is.na(bactgrowth)), aes(x = day, y = risk, group = bactgrowth)) +
  labs(y = 'Probability of colonisation by resistant Gram negative bacteria', 
       x = 'Antibiotic duration (days)') +
  facet_wrap(~ scenario) + 
  theme_minimal() + 
  theme(legend.position = 'bottom', 
        text = element_text(size = 20))

ggsave('~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/final_main/Rwithdays_treated.png',  
       dpi = 500, 
       width = 12, height = 7)


