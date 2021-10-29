rm(list = ls())


library(ggplot2); library(ggpubr)
library(boot)

get_metajags_plotdata <- function(saved_jagsoutput, days = 10){
  
  load(saved_jagsoutput)
  
  dat = jags.out[[1]]
  fit = jags.out[[2]]$BUGSoutput$summary
  arr = jags.out[[2]]$BUGSoutput$sims.matrix
  
  if (length(grep('norand', saved_jagsoutput)) > 0){
    
    ## get estimates for each data point 
    p = c()
    for (i in 1:max(lengths(dat))){ # each arm of any study 
      
      a = fit['a', '50%']
      b = fit['b', '50%']
      
      p[i] = inv.logit(a + #intercept
                         b * dat$abx_dur[i]) # slope for follow up period 
      
    }
    
    ## get predicted probabilities 
    plot.pred = data.frame(dur = 0:days, 
                           mean = NA,
                           low = NA,
                           high = NA)
    
    arr.q = t(apply(arr, 2, quantile, c(0.1, 0.25, 0.5, 0.75, 0.9)))
    
    for (i in 1:(days+1)){
      plot.pred$mean[i] = inv.logit(arr.q['a', '50%'] + arr.q['b', '50%'] * (i-1))
      
      plot.pred$lowlow[i] = inv.logit(arr.q['a', '10%'] + arr.q['b', '10%'] * (i-1))
      
      plot.pred$low[i] = inv.logit(arr.q['a', '25%'] + arr.q['b', '25%'] * (i-1))
      
      plot.pred$high[i] = inv.logit(arr.q['a', '75%'] + arr.q['b', '75%'] * (i-1))
      
      plot.pred$highhigh[i] = inv.logit(arr.q['a', '90%'] + arr.q['b', '90%'] * (i-1))
    }
    
    ## get data frame for plotting 
    plot.p = data.frame(dur = dat$abx_dur, 
                        y = p, 
                        actual = dat$n_ind_outcome,
                        sample_size = dat$n_ind_contributedsamples, 
                        study_id = as.factor(dat$pmid))
    
    RR.mean = exp(mean(arr.q['b', '50%']))
    RR10 = exp(mean(arr.q['b', '10%']))
    RR90 = exp(mean(arr.q['b', '90%']))
    
    
    
  } else {
    
    ## get estimates for each data point 
    p = c()
    for (i in 1:max(lengths(dat))){ # each arm of any study 
      
      study_id = dat$study_id[i]
      
      a = fit[paste0('a[', study_id, ']'), '50%']
      b = fit[paste0('b[', study_id, ']'), '50%']
      
      p[i] = inv.logit(a + #intercept
                         b * dat$abx_dur[i]) # slope for follow up period 
      
    }
    
    ## get predicted probabilities 
    plot.pred = data.frame(dur = 0:days, 
                           mean = NA,
                           low = NA,
                           high = NA)
    
    arr.q = t(apply(arr, 2, quantile, c(0.1, 0.25, 0.5, 0.75, 0.9)))
    
    for (i in 1:(days+1)){
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
    RR10 = exp(mean(arr.q[rownames(arr.q)[grep('b', rownames(arr.q))][1:dat$study_N], '10%']))
    RR90 = exp(mean(arr.q[rownames(arr.q)[grep('b', rownames(arr.q))][1:dat$study_N], '90%']))
    
  }
  
  RR = round(c(RR.mean, RR10, RR90), 2)
  
  return(list(plot.p, plot.pred, RR))
  
  
}

days = 15
dat = get_metajags_plotdata(saved_jagsoutput = 'runs/norand_0313nod_2021-10-27.Rdata', days = days)

mean = (dat[[3]][1])
low = (dat[[3]][2]) 
high = (dat[[3]][3]) 

dat[[1]]$actual = dat[[1]]$actual/dat[[1]]$sample_size

ggplot() + 
  geom_smooth(data = dat[[2]], aes(x = dur, y = mean), method = "lm", formula = y ~ splines::bs(x, 3), se = F, color = 'darkblue') + 
  geom_ribbon(data = dat[[2]], aes(x = dur, ymin = low, ymax = high), fill = 'grey', alpha = 0.4) + 
  geom_ribbon(data = dat[[2]], aes(x = dur, ymin = lowlow, ymax = highhigh), fill = 'grey', alpha = 0.2) + 
  geom_point(data = dat[[1]], aes(x = dur, y = actual, size = sample_size, color = study_id), alpha = 0.5) + 
  labs(x = 'Mean recorded antibiotic duration in each trial arm (days)', 
       y = 'Probability of resistant\nGram-negative bacteria colonisation') +
  # annotate('text', label = paste0('Increase in risk of acquiring resistance\ncarriage with one additional day of\nantibiotic is ', 
  #                                     mean, '% (80%CrI ', low, ' to ', high, '%)'), hjust = 0, x = 6.5, y = 0.1, size = 5) + 
  scale_y_continuous(limits = c(0 , 1), breaks = seq(0, 1, 0.1)) +
  #scale_x_continuous(limits = c(0 , days), breaks = seq(0, days, by = 2)) +
  scale_color_discrete(name = 'Study PMID') +
  scale_size_continuous(name = 'Number of participants') +
  theme_minimal() + 
  guides(colour = guide_legend(nrow = 2)) +
  theme(legend.position = 'bottom', 
        legend.box = "vertical",
        text = element_text(size = 18))

# odds ratio 
message(paste0('The odds ratio for being colonised with resistant bacteria per additional day of antibiotic treatment is ', mean, '(', low, ' to ', high, ').'))

ggsave('~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/final_main/meta.png',  
       dpi = 500,
       width = 12, height = 9)


############################################
#########SENSITIVITY ANALYSIS ##############
############################################
dat = get_metajags_plotdata(saved_jagsoutput = 'runs/randsens_0212nod_2021-10-13.Rdata', days = days)

ggplot() + 
  geom_smooth(data = dat[[2]], aes(x = dur, y = mean), method = "lm", formula = y ~ splines::bs(x, 3), se = F, color = 'darkblue') + 
  geom_ribbon(data = dat[[2]], aes(x = dur, ymin = low, ymax = high), fill = 'grey', alpha = 0.4) + 
  geom_ribbon(data = dat[[2]], aes(x = dur, ymin = lowlow, ymax = highhigh), fill = 'grey', alpha = 0.2) + 
  geom_point(data = dat[[1]], aes(x = dur, y = y, size = sample_size, color = study_id), alpha = 0.5) + 
  labs(x = 'Mean recorded antibiotic duration in each trial arm (days)', 
       y = 'Probability of resistant\nGram-negative bacteria colonisation') +
  # annotate('text', label = paste0('Increase in risk of acquiring resistance\ncarriage with one additional day of\nantibiotic is ', 
  #                                 mean, '% (80%CrI ', low, ' to ', high, '%)'), hjust = 0, x = 6.5, y = 0.1, size = 5) + 
  #scale_y_continuous(limits = c(0 , 1), breaks = seq(0, 1, 0.1)) +
  #scale_x_continuous(limits = c(0 , days), breaks = seq(0, days, by = 2)) +
  scale_color_discrete(name = 'Study PMID') +
  scale_size_continuous(name = 'Number of participants') +
  theme_minimal() + 
  guides(colour = guide_legend(nrow = 2)) +
  theme(legend.position = 'bottom', 
        legend.box = "vertical",
        text = element_text(size = 18))

