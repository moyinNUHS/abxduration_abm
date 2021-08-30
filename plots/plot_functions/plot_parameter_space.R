###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
###############################Check parameter space###############################
###################relationship between the parameters and outputs#################
###################################################################################

library(parallel)
library(ggplot2)
library(ggpubr)
library(reshape2)

source('msm_util_rtnorm.R')
source('los_abx_matrix.R')

### Ward charactersitics and antibiotic prescriptions
#####################################################

plot.los <- function (max.los.min, max.los.max, cum.r.1.min, cum.r.1.max, col) { #vary max.los, cum.r.1
  
  #number of runs to repeat 
  iterations = 3
  
  #all parameters that can affect los - to take the min and max of each parameter 
  n.bed = c(5, 50) 
  n.day = 300
  meanDur = c(3, 21)
  p.infect = c(0.1,1) 
  p.r.day1 = c(0.1,1)
  timestep = 1
  max.los = c(max.los.min, max.los.max)
  cum.r.1 = c(cum.r.1.min, cum.r.1.max)
  
  df = expand.grid(n.bed = n.bed, n.day = n.day, 
                   meanDur = meanDur, p.infect = p.infect, 
                   p.r.day1 = p.r.day1, timestep = timestep, 
                   max.los = max.los, cum.r.1 = cum.r.1)
  df = df[rep(row.names(df), each = iterations), ]
  df$runID = 1:nrow(df)
  
  feed.list = as.list(as.data.frame(t(df)))
  
  foo = function(x) {
    patient.matrix = los.abx.table(n.bed=x[1], n.day=x[2], meanDur=x[3], 
                                   p.infect=x[4], p.r.day1=x[5], timestep=x[6],
                                   max.los=x[7], cum.r.1=x[8])[[1]]
    los = table(patient.matrix)
  }
  
  save_runs = mclapply(feed.list, foo, mc.cores = 11)
  names(save_runs) = 1:nrow(df)
  plot.df = data.frame(runID = rep(names(save_runs), sapply(save_runs, length)),
                       los = unlist(save_runs))
  plot.df = merge(plot.df, df)
  
  if (col == "cum.r.1") {
    plot.df$group = ifelse(plot.df$cum.r.1 == cum.r.1.min, 'low', 'high')
  } else {
    plot.df$group = ifelse(plot.df$max.los == max.los.min, 'low', 'high')
  }
  
  p = ggplot(aes(x=los, group = runID, color = group),data = plot.df) + 
    geom_density(alpha=0.3) + 
    scale_colour_manual(values = c('orange', 'grey')) +
    ylab("Density") +
    xlab("Length of stay") + 
    annotate('text', label = paste('Orange represents', col, 'is at max'), x = 50, y = 0.1) +
    annotate('text', label = paste('Grey represents', col, 'is at min'), x = 50, y = 0.12) +
    labs(title=paste0('Effect of max.los and cum.r.1 on length of stay'), 
         subtitle = paste('p.infect and p.r.day1 0.1 : 1, meanDur 3:21, \nn.bed 5:50, max.los', max.los.min, ':', max.los.max, ", cum.r.1", cum.r.1.min, ":", cum.r.1.max)) +
    theme_minimal() +
    theme(legend.position = 'none') 
  
  return(p)
}

cum.r.1.plot <- function(cum.r.1.min, cum.r.1.max){
  
  #number of runs to repeat 
  iterations = 4
  
  #all parameters that can affect abx duration - to take the min and max of each parameter 
  n.bed = c(5, 50) 
  n.day = 300
  meanDur = c(3, 21)
  p.infect = 0 #no antibiotics started on admission
  p.r.day1 = 0 
  timestep = 1
  max.los = c(3, 30)
  cum.r.1 = seq(cum.r.1.min, cum.r.1.max, length.out = 10)
  
  df = expand.grid(n.bed = n.bed, n.day = n.day, 
                   meanDur = meanDur, p.infect = p.infect, 
                   p.r.day1 = p.r.day1, timestep = timestep, 
                   max.los = max.los, cum.r.1 = cum.r.1)
  df = df[rep(row.names(df), each = iterations), ]
  df$runID = 1:nrow(df)
  
  feed.list = as.list(as.data.frame(t(df)))
  
  foo = function(x) {
    abx.matrix = los.abx.table(n.bed=x[1], n.day=x[2], meanDur=x[3], 
                               p.infect=x[4], p.r.day1=x[5], timestep=x[6],
                               max.los=x[7], cum.r.1=x[8])[[2]]
    n.days.on.abx.r = sum(abx.matrix)/2/(dim(abx.matrix)[1] *dim(abx.matrix)[2])
  }
  
  save_runs = mclapply(feed.list, foo, mc.cores = 11)
  names(save_runs) = 1:nrow(df)
  plot.df = data.frame(runID = rep(names(save_runs), sapply(save_runs, length)),
                       prop.days.abxr = unlist(save_runs))
  plot.df = merge(df, plot.df)
  
  p = ggplot(aes(x= cum.r.1, y = prop.days.abxr), data = plot.df) + 
    geom_point(alpha=0.3, color = '#008080') + 
    ylab("Proportion of days on broad spectrum antibiotics\nin 300 days for ward size 5-50") +
    xlab("cum.r.1") + 
    labs(title=paste0('Effect of cum.r.1 on use of broad spectrum antibiotics'), 
         subtitle = paste('p.infect and p.r.day1 0, meanDur 3:21, nmax.los 3:30')) +
    theme_minimal() +
    theme(legend.position = 'none') 
  
  return(p)
  
}

### Baseline carriage status 
############################
plot.baseline <- function(model){
  
  #number of runs to repeat 
  iterations = 3
  intervals = 5
  
  #all parameters that can affect admission status - to take the min and max of each parameter 
  n.bed = 30 #does not matter 
  n.day = 300
  meanDur = c(3, 21)
  p.infect = 0.5 #does not matter 
  p.r.day1 = 0.5 #does not matter 
  timestep = 1
  max.los = 7 #does not matter
  cum.r.1 = 150 #does not matter 
  prop_R = seq(0.01, 0.8, length.out = intervals)
  
  if (model=='simple'){
    
    source('model_simple.R')
    prop_S = seq(0, 1, length.out = intervals)
    
    df = expand.grid(n.bed = n.bed, n.day = n.day, 
                     meanDur = meanDur, p.infect = p.infect, 
                     p.r.day1 = p.r.day1, timestep = timestep, 
                     max.los = max.los, cum.r.1 = cum.r.1, 
                     prop_S = prop_S, prop_R = prop_R)
    df = df[rep(row.names(df), each = iterations), ]
    
    feed.list = as.list(as.data.frame(t(df)))
    
    foo = function(x) {
      patient.matrix = los.abx.table(n.bed=x[1], n.day=x[2], meanDur=x[3], 
                                     p.infect=x[4], p.r.day1=x[5], timestep=x[6],
                                     max.los=x[7], cum.r.1=x[8])[[1]]
      los.array = summary.los(patient.matrix)
      admission = colo.table(patient.matrix, los.array, prop_R=x[10], prop_S=x[9])
      
      ss = length(which(admission == 'ss'))/max(patient.matrix)
      s = length(which(admission == 'S'))/max(patient.matrix)
      R = length(which(admission == 'R'))/max(patient.matrix)
      
      props.admission = c(ss=ss, s=s, R=R)
      
      return(props.admission)
    }
    
    save_runs = mclapply(feed.list, foo, mc.cores = 11)
    props = as.data.frame(matrix(unlist(save_runs), byrow = T, ncol = 3))
    colnames(props) = c('ss','S', 'R')
    plot.df = cbind(df, props)
    plot.dfm = melt(plot.df, measure.vars = c('ss', 'S', 'R'))
    
    pR = ggplot(aes(x = prop_R, y = value, color = variable), data = plot.dfm) + 
      geom_point(alpha = 0.3) +
      facet_grid(. ~ variable) +
      ylab('Proportion of carrier states on admission')+
      labs(title = 'Simple model', 
           subtitle = 'Effect of prop_R on admission carrier states') +
      theme_minimal()+
      theme(legend.position = 'none')
    
    pS = ggplot(aes(x = prop_S, y = value, color = variable), data = plot.dfm) + 
      geom_point(alpha = 0.3) +
      facet_grid(. ~ variable) +
      ylab('Proportion of carrier states on admission')+
      labs(title = 'Effect of prop_S on admission carrier states') +
      theme_minimal()+
      theme(legend.position = 'none')
   
    p = ggarrange(pR, pS, nrow = 2)
    
  } else if (model=='binary') {
    
    source('model_binary.R')
    
    prop_S = seq(0,1, length.out = intervals)
    prop_r = seq(0, 1, length.out = intervals)
    prop_Sr = seq(0,1, length.out = intervals)
    
    df = expand.grid(n.bed = n.bed, n.day = n.day, 
                     meanDur = meanDur, p.infect = p.infect, 
                     p.r.day1 = p.r.day1, timestep = timestep, 
                     max.los = max.los, cum.r.1 = cum.r.1, 
                     prop_S = prop_S, prop_R = prop_R, 
                     prop_r = prop_r, prop_Sr = prop_Sr)
    df = df[rep(row.names(df), each = iterations), ]
    
    feed.list = as.list(as.data.frame(t(df)))
    
    foo = function(x) {
      patient.matrix = los.abx.table(n.bed=x[1], n.day=x[2], meanDur=x[3], 
                                     p.infect=x[4], p.r.day1=x[5], timestep=x[6],
                                     max.los=x[7], cum.r.1=x[8])[[1]]
      los.array = summary.los(patient.matrix)
      admission = colo.table(patient.matrix, los.array, prop_R=x[10], prop_S=x[9], 
                             prop_r = x[11], prop_Sr = x[12])
      
      ss = length(which(admission == 'ss'))/max(patient.matrix)
      s = length(which(admission == 'S'))/max(patient.matrix)
      sR = length(which(admission == 'sR'))/max(patient.matrix)
      Sr = length(which(admission == 'Sr'))/max(patient.matrix)
      sr = length(which(admission == 'sr'))/max(patient.matrix)
      
      props.admission = c(ss=ss, s=s, sR=sR, Sr = Sr, sr = sr)
      
      return(props.admission)
    }
    
    save_runs = mclapply(feed.list, foo, mc.cores = 11)
    props = as.data.frame(matrix(unlist(save_runs), byrow = T, ncol = 5))
    colnames(props) = c('ss','S', 'sR', 'Sr' ,'sr')
    plot.df = cbind(df, props)
    plot.dfm = melt(plot.df, measure.vars = c('ss','S', 'sR', 'Sr' ,'sr'))
    
    pR = ggplot(aes(x = prop_R, y = value, color = variable), data = plot.dfm) + 
      geom_point(alpha = 0.3) +
      facet_grid(.~variable) +
      ylab('Proportion of carrier states on admission')+
      labs(title = 'Binary model', 
           subtitle = 'Effect of prop_R on admission carrier states') +
      theme_minimal()+
      theme(legend.position = 'bottom')
    
    pr = ggplot(aes(x = prop_r, y = value, color = variable), data = plot.dfm) + 
      geom_point(alpha = 0.3) +
      facet_grid(.~variable) +
      ylab('')+
      labs(title = 'Effect of prop_r on admission carrier states') +
      theme_minimal()+
      theme(legend.position = 'bottom')
    
    pS = ggplot(aes(x = prop_S, y = value, color = variable), data = plot.dfm) + 
      geom_point(alpha = 0.3) +
      facet_grid(.~variable) +
      ylab('')+
      labs(title = 'Effect of prop_S on admission carrier states') +
      theme_minimal()+
      theme(legend.position = 'bottom')
    
    pSr = ggplot(aes(x = prop_Sr, y = value, color = variable), data = plot.dfm) + 
      geom_point(alpha = 0.3) +
      facet_grid(.~variable) +
      ylab('Proportion of carrier states on admission')+
      labs(title = 'Effect of prop_Sr on admission carrier states') +
      theme_minimal()+
      theme(legend.position = 'bottom')
    
    p=ggarrange(pR, pr, pSr, pS, common.legend = T)

  } else if (model == 'frequency'){
    
    source('model_frequency.R')
    
    total_prop = c(0.1, 0.9)
    r_thres = c(6, 10)
    K = c(18, 24)
    
    df = expand.grid(n.bed = n.bed, n.day = n.day, 
                     meanDur = meanDur, p.infect = p.infect, 
                     p.r.day1 = p.r.day1, timestep = timestep, 
                     max.los = max.los, cum.r.1 = cum.r.1,  prop_R = prop_R, 
                     total_prop = total_prop, r_thres = r_thres, K = K)
    df = df[rep(row.names(df), each = iterations), ]
    
    feed.list = as.list(as.data.frame(t(df)))
    
    foo = function(x) {
      patient.matrix = los.abx.table(n.bed=x[1], n.day=x[2], meanDur=x[3], 
                                     p.infect=x[4], p.r.day1=x[5], timestep=x[6],
                                     max.los=x[7], cum.r.1=x[8])[[1]]
      los.array = summary.los(patient.matrix)
      admission = colo.table(patient.matrix, los.array, prop_R=x[9], total_prop=x[10],
                             r_thres=x[11], K = x[12])
      
      S = sum(exp(admission[[1]])/exp(admission[[3]]), na.rm = T)/max(patient.matrix)
      R = sum(exp(admission[[2]])/exp(admission[[3]]), na.rm = T)/max(patient.matrix)
      R.no = sum(admission[[2]] > x[11], na.rm = T)/max(patient.matrix)
      props.admission = c(S, R, R.no)
    
      return(props.admission)
    }
    
    save_runs = mclapply(feed.list, foo, mc.cores = 11)
    props = as.data.frame(matrix(unlist(save_runs), byrow = T, ncol = 3))
    colnames(props) = c('mean %S in total bact carried', 'mean %R in total bact carried', 'No. of R')
    plot.df = cbind(df, props)
    plot.dfm = melt(plot.df, measure.vars = c('mean %S in total bact carried', 'mean %R in total bact carried', 'No. of R'))
    
    p = ggplot(aes(x = prop_R, y = value, color = variable), data = plot.dfm) + 
      geom_point(alpha = 0.3) +
      facet_grid(.~variable) +
      ylab('')+
      labs(title = 'Frequency model', 
           subtitle = 'Effect of prop_R on admission carrier states') +
      theme_minimal()+
      theme(legend.position = 'none')
  }
  return (p)
}

### Within host dynamics 
########################
# model 2 repop.s and repop.r
plot.repop <- function (min, max, repop.name, intercept, los) {
  
  #what is the cumulative risk of repop with time 
  iter = 20
  repop = seq(min, max, length.out = iter)
  day.risk = data.frame(days = 1:los)
  
  for (i in 1:iter){
    daily.risk = repop[i] * (1 - repop[i])^(1:los-1) #if event did not happen in the previous days, and happens on a particular day 
    cum.daily.risk = cumsum(daily.risk)
    day.risk[,i+1] = cum.daily.risk
    colnames(day.risk)[i+1] = repop[i]
  }
  
  pd.melt = melt(day.risk, id.vars = 'days')
  
  p = ggplot(pd.melt, aes(x=days, y=value, group = variable, color=variable)) + 
    geom_line() + 
    scale_color_manual(values = c("red", rep("gray", iter-2), "blue")) +
    scale_fill_discrete(limits = c(min, max)) + 
    ylab('Cumulative daily probability of bacterial growth') +
    xlab('Length of stay (days)') + 
    theme_minimal() + 
    geom_vline(xintercept = intercept, linetype="dotted", color = "orange", size=1) +
    theme(legend.position = 'none') +
    annotate('text', x = 0, y = max(pd.melt$value)-0.1, label = paste0("Red represents ",repop.name, "=", min), hjust = 0) +
    annotate('text', x = 0, y = max(pd.melt$value), label = paste0("Blue represents ",repop.name, "=", max), hjust = 0)
  
  return(p)
}

#transmission pi_ssr 
#transmission, pi_ssr (models 1, 2, 3)
plot.pi_ssr <- function(min, max, los){
  
  ## for a certain los, what is the risk of transmission for various values of pi_ssr?
  iter=3
  pi_ssr = seq(min, max, length.out = iter)
  beds.r = seq(0.1, 1,length.out = iter)
  df = expand.grid(pi_ssr, beds.r)
  
  pd = data.frame(days = 1:los)
  for (i in 1:nrow(df)) {
    prop.transmit =  as.numeric(1-((1-df[i,][1])^(df[i,][2]))) #assume proportion of beds occupied by R is the same though out
    p.day= prop.transmit*(1-prop.transmit)^((1:los)-1) # get the probability for each day
    pd[, i+1] = cumsum(p.day)
    colnames(pd)[i+1] = paste0(i, df[i,][1]) #colnames of pd refers to pi_ssr values 
  } 
  
  pd.melt = melt(pd, id.vars = 'days')
  pd.melt$variable = as.character(pd.melt$variable)
  pd.melt$variable[grep(as.character(pi_ssr[1]),  pd.melt$variable)] =   pi_ssr[1]
  pd.melt$variable[grep(as.character(pi_ssr[2]),  pd.melt$variable)] =   pi_ssr[2]
  pd.melt$variable[grep(as.character(pi_ssr[3]),  pd.melt$variable)] =   pi_ssr[3]
  
  col1=alpha('red', 0.3)
  col2=alpha('gray', 0.3)
  col3=alpha('blue', 0.3)
  p = ggplot(pd.melt, aes(x=days, y=value, group = variable, color=variable)) + 
    geom_line() + 
    scale_color_manual(values = c(col1, col2, col3)) +
    scale_fill_discrete(limits = c(min, max)) + 
    ylab('Cumulative risk of transmission with duration of stay') +
    xlab('Length of stay (days)') + 
    theme_minimal() + 
    scale_x_continuous(breaks = seq(0,300, by =15))+
    theme(legend.position = 'none') +
    annotate('text', x = 100, y = 0.1, label = paste0("Red represents pi_ssr=", min, ", proportion of r bed=0.1"), hjust = 0) +
    annotate('text', x = 100, y = 0.16, label = paste0("Blue represents pi_ssr=", max, ", proportion of r bed=1"), hjust = 0)
  
  return (p)
}

# model 3 s_growth and r_growth 
growth <- function(min_r, max_r, min_s, max_s, los, paratoexplore){
  
  source('model_frequency.R')
  
  iter=2
  #min_r=0.01
  #max_r=0.05
  #min_s=0.005
  #max_s=0.015
  r_growth = seq(min_r, max_r, length.out = iter)
  s_growth = seq(min_s, max_s, length.out = iter)
  total_prop = seq(0.1, 0.9, length.out = iter)
  K = seq(18, 24, length.out = iter)
  r_thres = seq(6, 10, length.out = iter)
  prop_R = seq(0, 0.8, length.out = iter)
  abx.r = seq(1, 2, length.out = iter)
  df = expand.grid(r_growth=r_growth, total_prop=total_prop, prop_R= prop_R,K=K,r_thres=r_thres, s_growth=s_growth, abx.r = abx.r)
  
  n.bed = 30
  max.los = 300 #to make this same as the repop (above) to make comparison
  patient.matrix = matrix(rep(1:n.bed, each=max.los), ncol = n.bed)
  if (paratoexplore=='r_growth') {
    abx.matrix= matrix(rep(1, n.bed*max.los), ncol = n.bed)
  } else {
    abx.matrix= matrix(rep(0, n.bed*max.los), ncol = n.bed)
  }
  los.array=summary.los(patient.matrix)
  
  feed.list=list()
  for (i in 1:nrow(df)) {
    feed.list[[i]] = list(patient.matrix, los.array, abx.matrix, as.list(df[i,]))
  }
  
  foo = function(x) {
    colo.matrix = colo.table(patient.matrix=x[[1]], los.array=x[[2]], 
                             total_prop=as.numeric(x[[4]][2]), prop_R=as.numeric(x[[4]][3]), 
                             r_thres=as.numeric(x[[4]][5]), K = as.numeric(x[[4]][4]))
    filled.tab = nextDay(patient.matrix=x[[1]], los.array=x[[2]], abx.matrix=x[[3]], 
                         colo.matrix = colo.matrix, pi_ssr=0, total_prop=as.numeric(x[[4]][2]),  
                         K = as.numeric(x[[4]][4]), r_growth = as.numeric(x[[4]][1]), r_thres=as.numeric(x[[4]][5]), 
                         r_trans=0, s_growth = as.numeric(x[[4]][6]), abx.s=0, abx.r=as.numeric(x[[4]][7]), timestep=1)
    list(grow_s = exp(filled.tab[[1]]), grow_r = exp(filled.tab[[2]]))
  }
  
  save_runs = mclapply(feed.list, foo, mc.cores = 11)
  
  s.p.list=r.p.list=list()
  for (i in 1:nrow(df)){
    s.df = as.data.frame(save_runs[[i]][1])
    s.df$days = 1:300
    colnames(s.df)[1:n.bed] = paste0(i,colnames(s.df)[1:n.bed])
    s.df.m=melt(s.df, id.vars = 'days')
    s.df.m$r_growth = rep(unlist(df[i,][1]), nrow(s.df.m))
    s.df.m$total_prop = rep(unlist(df[i,][2]), nrow(s.df.m))
    s.df.m$prop_R = rep(unlist(df[i,][3]), nrow(s.df.m))
    s.df.m$K = rep(unlist(df[i,][4]), nrow(s.df.m))
    s.df.m$r_thres = rep(unlist(df[i,][5]), nrow(s.df.m))
    s.df.m$s_growth = rep(unlist(df[i,][6]), nrow(s.df.m))
    s.df.m$abx.r = rep(unlist(df[i,][7]), nrow(s.df.m))
    
    r.df = as.data.frame(save_runs[[i]][2])
    r.df$days = 1:300
    colnames(r.df)[1:n.bed] = paste0(i,colnames(r.df)[1:n.bed])
    r.df.m=melt(r.df, id.vars = 'days')
    r.df.m$r_growth = rep(unlist(df[i,][1]), nrow(r.df.m))
    r.df.m$total_prop = rep(unlist(df[i,][2]), nrow(r.df.m))
    r.df.m$prop_R = rep(unlist(df[i,][3]), nrow(r.df.m))
    r.df.m$K = rep(unlist(df[i,][4]), nrow(r.df.m))
    r.df.m$r_thres = rep(unlist(df[i,][5]), nrow(r.df.m))
    r.df.m$s_growth = rep(unlist(df[i,][6]), nrow(r.df.m))
    r.df.m$abx.r = rep(unlist(df[i,][7]), nrow(r.df.m))
    
    r.p.list[[i]]= r.df.m
    s.p.list[[i]]= s.df.m
  }
  
  r.pd = do.call('rbind', r.p.list)
  s.pd = do.call('rbind', s.p.list)
  
  colnumber= which(colnames(s.pd)==paratoexplore)
  min.value= unique(s.pd[colnumber])[1,] #min value of the parameter to explore
  max.value= unique(s.pd[colnumber])[2,] #max value of the parameter to explore
  pd.plot1s = s.pd[which(s.pd[colnumber]==min.value),] #df with min value of the parameter
  pd.plot1r = r.pd[which(r.pd[colnumber]==min.value),]
  pd.plot2s = s.pd[which(s.pd[colnumber]==max.value),]
  pd.plot2r = r.pd[which(r.pd[colnumber]==max.value),]
  min.r = paste0(round(min (r.pd$value)/1000000000), '(10^9)')
  max.r = paste0(round(max (r.pd$value)/1000000000), '(10^9)')
  min.s = paste0(round(min (s.pd$value)/1000000000), '(10^9)')
  max.s = paste0(round(max (s.pd$value)/1000000000), '(10^9)')
  
  col=alpha('red', 0.4)
  p1s = ggplot(pd.plot1s, aes(x=days, y=value, color=variable)) + 
    geom_line() + ylab('No. of S bacteria') + 
    scale_color_manual(values = rep(col, (nrow(pd.plot1s)/300)), guide =F) +
    xlab('Length of stay (days)') + labs(title=paste0('No. of bacteria with varying prop_R, K, r_growth, s_growth, total_prop, r_thres\n',
                                                      'min of R=', min.r, ' max of R=', max.r,' min of S=', min.s,' max of S=', max.s),
                                         subtitle= paste0(paratoexplore, '=', min.value)) + theme_minimal() +
    theme(legend.position = 'bottom')
  p2s = ggplot(pd.plot2s, aes(x=days, y=value, color=variable)) + 
    geom_line() + ylab('No. of S') +
    scale_color_manual(values = rep(col, (nrow(pd.plot2s)/300)), guide =F) +
    xlab('Length of stay (days)') + labs(title='  \n  ',
                                         subtitle= paste0(paratoexplore, '=', max.value)) + theme_minimal()+
    theme(legend.position = 'bottom')
  p1r = ggplot(pd.plot1r, aes(x=days, y=value, color=variable)) + 
    geom_line() + ylab('No. of R bacteria') + theme_minimal()+
    scale_color_manual(values = rep(col, (nrow(pd.plot2s)/300)), guide =F) +
    xlab('Length of stay (days)') + theme(legend.position = 'bottom')
  p2r = ggplot(pd.plot2r, aes(x=days, y=value, color=variable)) + 
    geom_line() + ylab('No. of R') + theme_minimal()+
    scale_color_manual(values = rep(col, (nrow(pd.plot2s)/300)), guide =F) +
    xlab('Length of stay (days)') + theme(legend.position = 'bottom')
  
  p=ggarrange(p1s, p2s, p1r, p2r, ncol=2, nrow = 2, common.legend = T)
  
  return (p)
}

#antibiotic killing 
# model 1 and 2 
abx.kill12 <- function (min, max, los, abx.name) {
  
  #what is the cumulative risk of abx killing with increasing duration of stay? 
  iter=20
  abx=seq(min, max, length.out = iter)
  day.risk=data.frame(days=1:los)
  
  for (i in 1:iter){
    daily.risk = abx[i] * (1 - abx[i])^(1:los-1) #if event did not happen in the previous days, and happened on a particular day 
    cum.daily.risk = cumsum(daily.risk)
    day.risk[,i+1] = cum.daily.risk
    colnames(day.risk)[i+1] = abx[i]
  }
  
  pd.melt = melt(day.risk, id.vars = 'days')
  
  p= ggplot(pd.melt, aes(x=days, y=value, group = variable, color=variable)) + 
    geom_line() + 
    scale_color_manual(values = c("red", rep("gray", iter-2), "blue")) +
    scale_fill_discrete(limits = c(min, max)) + 
    ylab('Cumulative risk of bactrial killing with duration of stay') +
    xlab('Length of stay (days)') + 
    theme_minimal() + 
    theme(legend.position = 'none') +
    annotate('text', x = 3, y = min(pd.melt$value), label = paste0("Red represents ",abx.name, "=", min), hjust = 0) +
    annotate('text', x = 3, y = min(pd.melt$value)+0.1, label = paste0("Blue represents ",abx.name, "=", max), hjust = 0)
  
  return(p)
}

# model 3 
abx.kill3 <- function(min_r, max_r, min_s, max_s, los, abx.name, abx_min, abx_max, paratoexplore){
  
  source('model_frequency.R')
  iter=2
  # min_r=0.01
  # max_r=0.05
  # min_s=0.005
  # max_s=0.015
  r_growth = seq(min_r, max_r, length.out = iter)
  s_growth = seq(min_s, max_s, length.out = iter)
  total_prop = seq(0.1, 0.9, length.out = iter)
  K = seq(18, 24, length.out = iter)
  r_thres = seq(6, 10, length.out = iter)
  prop_R = seq(0, 0.8, length.out = iter)
  abx = seq (abx_min, abx_max, length.out = iter)
  df = expand.grid(r_growth=r_growth, total_prop=total_prop, prop_R= prop_R,K=K,r_thres=r_thres, s_growth=s_growth, abx=abx)
  
  n.bed = 30
  max.los = los 
  patient.matrix = matrix(rep(1:n.bed, each=max.los), ncol = n.bed)
  if (abx.name == 'abx.r') {
    abx.matrix= matrix(rep(2, n.bed*max.los), ncol = n.bed)
  } else if (abx.name == 'abx.s') {
    abx.matrix= matrix(rep(1, n.bed*max.los), ncol = n.bed)
  }
  los.array=summary.los(patient.matrix)
  
  feed.list=list()
  for (i in 1:nrow(df)) {
    feed.list[[i]] = list(patient.matrix, los.array, abx.matrix, as.list(df[i,]))
  }
  
  foo = function(x) {
    colo.matrix = colo.table(patient.matrix=x[[1]], los.array=x[[2]], 
                             total_prop=as.numeric(x[[4]][2]), prop_R=as.numeric(x[[4]][3]), 
                             r_thres=as.numeric(x[[4]][5]), K = as.numeric(x[[4]][4]))
    filled.tab = nextDay(patient.matrix=x[[1]], los.array=x[[2]], abx.matrix=x[[3]], 
                         colo.matrix = colo.matrix, pi_ssr=0, total_prop=as.numeric(x[[4]][2]),  
                         K = as.numeric(x[[4]][4]), r_growth = as.numeric(x[[4]][1]), r_thres=as.numeric(x[[4]][5]), 
                         r_trans=0, s_growth = as.numeric(x[[4]][6]), abx.s=as.numeric(x[[4]][7]), abx.r=as.numeric(x[[4]][7]), timestep=1)
    list(grow_s = exp(filled.tab[[1]]), grow_r = exp(filled.tab[[2]]))
  }
  
  save_runs = mclapply(feed.list, foo, mc.cores = 11)
  
  s.p.list=r.p.list=list()
  for (i in 1:nrow(df)){
    s.df = as.data.frame(save_runs[[i]][1])
    s.df$days = 1:los
    colnames(s.df)[1:n.bed] = paste0(i,colnames(s.df)[1:n.bed])
    s.df.m=melt(s.df, id.vars = 'days')
    s.df.m$r_growth = rep(unlist(df[i,][1]), nrow(s.df.m))
    s.df.m$total_prop = rep(unlist(df[i,][2]), nrow(s.df.m))
    s.df.m$prop_R = rep(unlist(df[i,][3]), nrow(s.df.m))
    s.df.m$K = rep(unlist(df[i,][4]), nrow(s.df.m))
    s.df.m$r_thres = rep(unlist(df[i,][5]), nrow(s.df.m))
    s.df.m$s_growth = rep(unlist(df[i,][6]), nrow(s.df.m))
    s.df.m$abx= rep(unlist(df[i,][7]), nrow(s.df.m))
    
    r.df = as.data.frame(save_runs[[i]][2])
    r.df$days = 1:los
    colnames(r.df)[1:n.bed] = paste0(i,colnames(r.df)[1:n.bed])
    r.df.m=melt(r.df, id.vars = 'days')
    r.df.m$r_growth = rep(unlist(df[i,][1]), nrow(r.df.m))
    r.df.m$total_prop = rep(unlist(df[i,][2]), nrow(r.df.m))
    r.df.m$prop_R = rep(unlist(df[i,][3]), nrow(r.df.m))
    r.df.m$K = rep(unlist(df[i,][4]), nrow(r.df.m))
    r.df.m$r_thres = rep(unlist(df[i,][5]), nrow(r.df.m))
    r.df.m$s_growth = rep(unlist(df[i,][6]), nrow(r.df.m))
    r.df.m$abx = rep(unlist(df[i,][7]), nrow(r.df.m))
    
    r.p.list[[i]]= r.df.m
    s.p.list[[i]]= s.df.m
  }
  
  r.pd = do.call('rbind.data.frame', r.p.list)
  s.pd = do.call('rbind.data.frame', s.p.list)
  
  colnumber= which(colnames(s.pd)==paratoexplore)
  min.value= unique(s.pd[colnumber])[1,]
  max.value= unique(s.pd[colnumber])[2,]
  pd.plot1s = s.pd[which(s.pd[colnumber]==min.value),]
  pd.plot1r = r.pd[which(r.pd[colnumber]==min.value),]
  pd.plot2s = s.pd[which(s.pd[colnumber]==max.value),]
  pd.plot2r = r.pd[which(r.pd[colnumber]==max.value),]
  min.r = paste0(round(min (r.pd$value)/1000000000), '(10^9)')
  max.r = paste0(round(max (r.pd$value)/1000000000), '(10^9)')
  min.s = paste0(round(min (s.pd$value)/1000000000), '(10^9)')
  max.s = paste0(round(max (s.pd$value)/1000000000), '(10^9)')
  
  col=alpha('red', 0.4)
  p1s = ggplot(pd.plot1s, aes(x=days, y=value, color=variable)) + 
    geom_line() + ylab('No. of S bacteria') + 
    scale_x_continuous(breaks=seq(0,los,by=3))+
    scale_color_manual(values = rep(col, (nrow(pd.plot1s)/los)), guide =F) +
    xlab('Length of stay (days)') + labs(title=paste0('No. of bacteria with varying prop_R, K, r_growth, s_growth, total_prop, r_thres\n',
                                                      'min of R=', min.r, ' max of R=', max.r,' min of S=', min.s,' max of S=', max.s),
                                         subtitle= paste0(paratoexplore, '=', min.value)) + theme_minimal() +
    theme(legend.position = 'bottom')
  p2s = ggplot(pd.plot2s, aes(x=days, y=value, color=variable)) + 
    geom_line() + ylab('No. of S') +
    scale_x_continuous(breaks=seq(0,los,by=3))+
    scale_color_manual(values = rep(col, (nrow(pd.plot2s)/los)), guide =F) +
    xlab('Length of stay (days)') + labs(title='  \n  ',
                                         subtitle= paste0(paratoexplore, '=', max.value)) + theme_minimal()+
    theme(legend.position = 'bottom')
  p1r = ggplot(pd.plot1r, aes(x=days, y=value, color=variable)) + 
    geom_line() + ylab('No. of R bacteria') + theme_minimal()+
    scale_x_continuous(breaks=seq(0,los,by=3))+
    scale_color_manual(values = rep(col, (nrow(pd.plot2s)/los)), guide =F) +
    xlab('Length of stay (days)') + theme(legend.position = 'bottom')
  p2r = ggplot(pd.plot2r, aes(x=days, y=value, color=variable)) + 
    geom_line() + ylab('No. of R') + theme_minimal()+
    scale_x_continuous(breaks=seq(0,los,by=3))+
    scale_color_manual(values = rep(col, (nrow(pd.plot2s)/los)), guide =F) +
    xlab('Length of stay (days)') + theme(legend.position = 'bottom')
  
  p=ggarrange(p1s, p2s, p1r, p2r, ncol=2, nrow = 2, common.legend = T)
  
  return (p)
}

###decolonisation 
##################################
#decolonisation, mu (models 1 and 2) 
plot.mu <- function(min, max, los){
  
  ## for a certain los, what is the risk of decolonisation for various values of mu?
  iter=20
  mu = seq(min, max, length.out = iter)
  
  pd = data.frame(days = 1:los)
  for (i in 1:iter) {
    p.day= mu[i]*(1-mu[i])^((1:los)-1) # get the probability for each day
    pd[, i+1] = cumsum(p.day)
    colnames(pd)[i+1] = mu[i]
  }
  
  pd.melt = melt(pd, id.vars = 'days')
  
  p = ggplot(pd.melt, aes(x=days, y=value, group = variable, color=variable)) + 
    geom_line() + 
    scale_color_manual(values = c("red", rep("gray", iter-2), "blue")) +
    scale_fill_discrete(limits = c(min, max)) + 
    ylab('Cumulative risk of decolonisation with duration of stay') +
    xlab('Length of stay (days)') + 
    theme_minimal() + 
    theme(legend.position = 'none') +
    annotate('text', x = 50, y = 0.2, label = paste0("Red represents mu=", min), hjust = 0) +
    annotate('text', x = 50, y = 0.26, label = paste0("Blue represents mu=", max), hjust = 0)
  
  return (p)
}



