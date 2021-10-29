dat.to.analyse1 = dat.to.analyse

dat.to.analyse = dat.to.analyse[order(dat.to.analyse$abx_dur),]
dat.to.analyse = dat.to.analyse[order(dat.to.analyse$PMID),]

mean = dat.to.analyse$n_ind_outcome/ dat.to.analyse$n_ind_contributedsamples

lower = upper = c()
for (i in 1:nrow(dat.to.analyse)){
  pt = prop.test(dat.to.analyse$n_ind_outcome[i], dat.to.analyse$n_ind_contributedsamples[i])
  lower[i] = pt$conf.int[1]
  upper[i] =  pt$conf.int[2]
}

d = data.frame(mean, lower, upper)
d = rbind(rep(NA, 3), d)
tabletext = dat.to.analyse[,c('PMID','bacteria_reported', 'abx_dur', 'resistance_type', 'time_baseline_endoffu')]
tabletext$PMID[duplicated(tabletext$PMID)] = ''
tabletext = rbind(colnames(tabletext), tabletext)

d %>% 
  forestplot(labeltext = tabletext, 
             is.summary = rep(FALSE, nrow(d)+1),
             clip = c(0.1, 2.5), 
             graph.pos=6,
             graphwidth=unit (10, "cm"),
             xlog = F, 
             col = fpColors(box = "royalblue",
                            line = "darkblue",
                            summary = "royalblue"))
