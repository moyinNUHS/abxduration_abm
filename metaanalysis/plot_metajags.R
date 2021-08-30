rm(list = ls())


library(ggplot2); library(ggpubr)
library(boot)

get_metajags_plotdata <- function(saved_jagsoutput){
  
  load(saved_jagsoutput)
  
  dat = jags.out[[1]]
  fit = jags.out[[2]]$BUGSoutput$summary
  arr = jags.out[[2]]$BUGSoutput$sims.matrix
  
  ## get estimates for each data point 
  p = c()
  for (i in 1:max(lengths(dat))){
    
    study_id = dat$study_id[i]
    
    a = fit[paste0('a[', study_id, ']'), '50%']
    b = fit[paste0('b[', study_id, ']'), '50%']
    c = ifelse(i %in% dat$binom_dist, fit[paste0('c[', study_id, ']'), '50%'], 0)
    
    p[i] = inv.logit(a + #intercept
                       b * dat$abx_dur[i] + # slope for duration
                       c * dat$time_baseline_endoffu[i] ) # slope for follow up period 
    
    
  }
  
  ## get predicted probabilities 
  plot.pred = data.frame(dur = 0:20, 
                    mean = NA,
                    low = NA,
                    high = NA)
  
  arr.q = t(apply(arr, 2, quantile, c(0.1, 0.25, 0.5, 0.75, 0.9)))
  
  for (i in 1:21){
    plot.pred$mean[i] = inv.logit(mean(fit[rownames(fit)[grep('a', rownames(fit))][1:dat$study_N], '50%']) + 
                          mean(fit[rownames(fit)[grep('b', rownames(fit))][1:dat$study_N], '50%']) * (i-1))
    
    plot.pred$lowlow[i] = inv.logit(mean(arr.q[rownames(arr.q)[grep('a', rownames(arr.q))][1:dat$study_N], '10%']) + 
                                   mean(arr.q[rownames(arr.q)[grep('b', rownames(arr.q))][1:dat$study_N], '10%'])* (i-1))
    
    plot.pred$low[i] = inv.logit(mean(arr.q[rownames(arr.q)[grep('a', rownames(arr.q))][1:dat$study_N], '25%']) + 
                         mean(arr.q[rownames(arr.q)[grep('b', rownames(arr.q))][1:dat$study_N], '25%'])* (i-1))
    
    plot.pred$high[i] = inv.logit(mean(arr.q[rownames(arr.q)[grep('a', rownames(arr.q))][1:dat$study_N], '75%']) + 
                         mean(arr.q[rownames(arr.q)[grep('b', rownames(arr.q))][1:dat$study_N], '75%'])* (i-1))
    
    plot.pred$highhigh[i] = inv.logit(mean(arr.q[rownames(arr.q)[grep('a', rownames(arr.q))][1:dat$study_N], '90%']) + 
                                    mean(arr.q[rownames(arr.q)[grep('b', rownames(arr.q))][1:dat$study_N], '90%'])* (i-1))
  }
  
  ## get data frame for plotting 
  plot.p = data.frame(dur = dat$abx_dur, 
                      y = p, 
                      sample_size = dat$n_ind_contributedsamples, 
                      study_id = as.factor(dat$pmid))
  
  RR.mean = exp(mean(arr.q[rownames(arr.q)[grep('b', rownames(arr.q))][1:dat$study_N], '50%']))
  RR.2.5 = exp(mean(arr.q[rownames(arr.q)[grep('b', rownames(arr.q))][1:dat$study_N], '10%']))
  RR.97.5 = exp(mean(arr.q[rownames(arr.q)[grep('b', rownames(arr.q))][1:dat$study_N], '90%']))
  
  RR = round(c(RR.mean, RR.2.5, RR.97.5), 2)
  
  return(list(plot.p, plot.pred, RR))
  
}

np.dat = get_metajags_plotdata(saved_jagsoutput = 'metaanalysis/runs/jagsOUT_noappro_30d_all_2021-06-28.Rdata')
p.dat = get_metajags_plotdata(saved_jagsoutput = 'metaanalysis/runs/jagsOUT_appro_30d_all_2021-06-28.Rdata')

plot.d = rbind(np.dat[[1]], p.dat[[1]])
plot.d$type = rep(c('Resistance phenotypes treated with ineffective antibiotic', 'Resistance phenotypes treated with effective antibiotic'), c(nrow(np.dat[[1]]), nrow(p.dat[[1]])))
plot.d$type = as.factor(plot.d$type)
plot.d$type = factor(plot.d$type, levels = c('Resistance phenotypes treated with ineffective antibiotic', 'Resistance phenotypes treated with effective antibiotic'))


plot.pred = rbind(np.dat[[2]], p.dat[[2]])
plot.pred$type = rep(c('Resistance phenotypes treated with ineffective antibiotic', 'Resistance phenotypes treated with effective antibiotic'), c(nrow(np.dat[[2]]), nrow(p.dat[[2]])))
plot.pred$type = as.factor(plot.pred$type)
plot.pred$type = factor(plot.pred$type, levels = c('Resistance phenotypes treated with ineffective antibiotic', 'Resistance phenotypes treated with effective antibiotic'))

rr.text = data.frame(dur = 5, y = 0.8,lab = "Text",
                     mean = c(np.dat[[3]][1], p.dat[[3]][1]),
                     low = c(np.dat[[3]][2], p.dat[[3]][2]),
                     high = c(np.dat[[3]][3], p.dat[[3]][3])
)
rr.text$type = c('Resistance phenotypes treated with ineffective antibiotic', 'Resistance phenotypes treated with effective antibiotic')
rr.text$type = as.factor(rr.text$type)
rr.text$type = factor(rr.text$type, levels = c('Resistance phenotypes treated with ineffective antibiotic', 'Resistance phenotypes treated with effective antibiotic'))


ggplot() + 
  geom_point(data = plot.d, aes(x = dur, y = y, size = sample_size, color = study_id), alpha = 0.5) + 
  geom_smooth(data = plot.pred, aes(x = dur, y = mean, group = type), method = "lm", formula = y ~ splines::bs(x, 15), se = F, color = 'darkblue') + 
  geom_ribbon(data = plot.pred, aes(x = dur, ymin = low, ymax = high, group = type), fill = 'grey', alpha = 0.4) + 
  geom_ribbon(data = plot.pred, aes(x = dur, ymin = lowlow, ymax = highhigh, group = type), fill = 'grey', alpha = 0.2) + 
  labs(x = 'Antibiotic duration (days)', y = 'Probability of colonisation by resistant Gram negative bacteria') +
  geom_text(data = rr.text, aes(x = dur, y = y, label = paste0('Daily increase in risk of acquiring resistance\ncarriage with one additional day of antibiotic\nis ', (mean-1)*100, '% (80%CrI ', (low-1)*100, ' to ', (high-1)*100, '%)')), hjust = 0) + 
  scale_y_continuous(limits = c(0 , 1), breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(limits = c(0 , 20)) +
  scale_color_discrete(name = 'Study PMID') +
  scale_size_continuous(name = 'Number of participants') +
  facet_wrap(~ type) +
  theme_minimal() + 
  guides(colour = guide_legend(nrow = 3)) +
  theme(legend.position = 'bottom', 
        legend.box = "vertical",
        text = element_text(size = 18))

ggsave('~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/final_main/meta.png',  
       dpi = 500,
       width = 12, height = 9)

####### 
#sensitivity analyses

(np.60 = get_metajags_plotdata(saved_jagsoutput = 'runs/jagsOUT_noappro_60d_all_2021-06-28.Rdata')[[3]])
(np.hosp = get_metajags_plotdata(saved_jagsoutput = 'runs/jagsOUT_noappro_30d_hosp_2021-06-28.Rdata')[[3]])

(p.60 = get_metajags_plotdata(saved_jagsoutput = 'runs/jagsOUT_appro_60d_all_2021-06-28.Rdata')[[3]])
(p.hosp = get_metajags_plotdata(saved_jagsoutput = 'runs/jagsOUT_appro_30d_hosp_2021-06-28.Rdata')[[3]])


