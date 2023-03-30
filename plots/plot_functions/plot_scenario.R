#########################################################################
#######Effect of antibiotic duration on hospitalised patients############
#####################Show one run per scenario ##########################
#########################################################################

#### CLEAN DATA FOR PLOTTING 

dataformosaic <- function(data, label, n.bed = para$n.bed, n.day = para$n.day, timestep=1){
    
    data.t=t(data)
    rownames(data.t)=as.factor(1:n.bed)
    colnames(data.t)=as.factor(1:(n.day*timestep))
    data.melt=melt(data.t)
    colnames(data.melt)=c('Bed number','Time (days)',label)
    
    return(data.melt)
    
}

admitdays <- function(patient.matrix){
    
    change = apply(patient.matrix, 2, function (v) c(1, 1 + which(diff(v) != 0)) - 0.5)
    x = x.end = as.vector(unlist(change))
    y = rep(1:ncol(patient.matrix), lengths(change))-0.5
    df = data.frame(x=x, x.end=x.end, y=y, y.end=y+1)
    
    return(df)
}

#### RUN MODEL 

run_scenario <- function(para, n.runs) {
    
    #Run model 
    sdDur = 1
    
    store = list()
    
    for (i in 1:n.runs){
        matrixes = los_abx_table_varydur(n.bed = para$n.bed, n.day = para$n.day, max.los = para$max.los, 
                                         p.infect = para$p.infect, p.r.day1 = para$p.r.day1, p.r.after = para$p.r.after, p.infect.after = para$p.infect.after, 
                                         meanDur = para$short_dur, timestep = 1)
        patient.matrix.short = matrixes[[1]]
        day1.short = admitdays(patient.matrix.short)
        
        abx.matrix.short = matrixes[[2]]
        los.array.short = summary_los(patient.matrix=patient.matrix.short)
        colo.matrix.short = colo.table(patient.matrix=patient.matrix.short, los=los.array.short, 
                                       prop_R=para$prop_R, prop_r=para$prop_r, prop_Sr=para$prop_Sr, prop_S=para$prop_S)
        
        colo_table_filled_short = nextDay(patient.matrix=patient.matrix.short, abx.matrix=abx.matrix.short, colo.matrix=colo.matrix.short, 
                                          pi_ssr=para$pi_ssr, bif=para$bif, mu=para$mu, fitness.r=para$fitness.r,
                                          repop.s=para$repop.s, abx.r=para$abx.r, abx.s=para$abx.s, timestep=1)
        
        matrixes = los_abx_table_varydur(n.bed=para$n.bed, n.day=para$n.day, max.los=para$max.los, 
                                         p.infect=para$p.infect, p.r.day1=para$p.r.day1, p.infect.after=para$p.infect.after, 
                                         p.r.after = para$p.r.after,
                                         meanDur= para$long_dur, timestep=1)
        patient.matrix.long=matrixes[[1]]
        day1.long = admitdays(patient.matrix.long)
        
        abx.matrix.long=matrixes[[2]]
        los.array.long = summary_los(patient.matrix=patient.matrix.long)
        colo.matrix.long = colo.table(patient.matrix=patient.matrix.long, los=los.array.long, 
                                      prop_R=para$prop_R, prop_r=para$prop_r, prop_Sr=para$prop_Sr, prop_S=para$prop_S)
        
        colo_table_filled_long = nextDay(patient.matrix=patient.matrix.long, abx.matrix=abx.matrix.long, colo.matrix=colo.matrix.long, 
                                         pi_ssr=para$pi_ssr, bif=para$bif, mu=para$mu, fitness.r=para$fitness.r,
                                         repop.s=para$repop.s, abx.r=para$abx.r, abx.s=para$abx.s, timestep=1)
        
        if (i == 1) {
            store[['abx.matrix.short']] = abx.matrix.short
            store[['abx.matrix.long']] = abx.matrix.long
            store[['day1.long']] = day1.long
            store[['day1.short']] = day1.short
        }
        
        store[['colo_table_filled_short']][[i]] = colo_table_filled_short
        store[['colo_table_filled_long']][[i]] = colo_table_filled_long
    }
    
    ## to produce plots need 
    # Plot 1: 1 * abx.matrix.short, abx.matrix.long
    # Plot 2: 1 * colo_table_filled_short, colo_table_filled_long 
    # Plot 3 and 4: multiple * colo_table_filled_short, colo_table_filled_long 
    
    return(store)
    
}

##### PLOT 

winteralmond = c('white',"#87C2BE","#5E8E7B")
lgdcol = rgb(0.185, 0.188, 0.154, alpha = .05)
base_size = 9

plot_scenario <- function(para, dat) {
    
    ### PLOT 1 ###
    ####Abx use plots 
    abx.matrix.short = dat$abx.matrix.short
    day1.short = dat$day1.short
    mosaicdata.abx.short = dataformosaic(data=abx.matrix.short ,label='Antibiotic type',n.bed=para$n.bed, n.day=para$n.day, timestep=1)
    mosaicdata.abx.short$`Antibiotic type`=as.factor(mosaicdata.abx.short$`Antibiotic type`)
    abx.matrix.long = dat$abx.matrix.long
    day1.long = dat$day1.long
    mosaicdata.abx.long = dataformosaic(data=abx.matrix.long ,label='Antibiotic type',n.bed=para$n.bed, n.day=para$n.day, timestep=1)
    mosaicdata.abx.long$`Antibiotic type`=as.factor(mosaicdata.abx.long$`Antibiotic type`)
    
    p.abx.short = ggplot(mosaicdata.abx.short, aes(`Time (days)`, `Bed number`)) + 
        geom_tile(aes(fill = `Antibiotic type`)) + 
        scale_fill_manual(values=winteralmond,labels = c("None", "Effective only against susceptible organisms", "Effective against susceptible and resistant organisms"), name='Antibiotic type')+
        theme_bw(base_size=base_size) + 
        labs(title= 'Short antibiotic treatment duration',
             x = "Observation days", y = "Unique patients")+ 
        scale_x_continuous(expand = c(0, 0), breaks= seq(10, para$n.day, 10), labels = as.character(seq(10, para$n.day, 10))) +
        scale_y_continuous(expand = c(0, 0), breaks= 1:10) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.title = element_text(size = base_size+2), 
              legend.text  = element_text(size = base_size+2),
              axis.text = element_text(size = base_size+2, colour = "grey50"), 
              axis.title = element_text(size=base_size+2))+
        geom_segment(data=day1.short, aes(x,y,xend=x.end, yend=y.end), size=0.25, inherit.aes=F, colour='blue')
    
    p.abx.long = ggplot(mosaicdata.abx.long, aes(`Time (days)`, `Bed number`)) + 
        geom_tile(aes(fill = `Antibiotic type`)) + 
        scale_fill_manual(values=winteralmond,labels = c("None", "Effective only against susceptible organisms", "Effective against susceptible and resistant organisms"), name='Antibiotic type')+
        theme_bw(base_size=base_size) + 
        labs(title= 'Long antibiotic treatment duration',
             x = "Observation days", y="")+ 
        scale_x_continuous(expand = c(0, 0), breaks= seq(10, para$n.day, 10), labels = as.character(seq(10, para$n.day, 10))) +
        scale_y_continuous(expand = c(0, 0), breaks= 1:10) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.title = element_text(size = base_size+2), 
              legend.text  = element_text(size = base_size+2),
              axis.text = element_text(size = base_size+2, colour = "grey50"), 
              axis.title = element_text(size=base_size+2))+
        geom_segment(data=day1.long, aes(x,y,xend=x.end, yend=y.end), size=0.25, inherit.aes=F, colour='blue')
    
    
    abx.mosaic=ggarrange(p.abx.short, p.abx.long, ncol=2, common.legend = T, legend = 'bottom')
    
    ### PLOT 2 ###
    ##Carriage mosaic 
    colo_table_filled_short = dat$colo_table_filled_short[[1]]
    mosaicdata.car.short = dataformosaic(data=colo_table_filled_short,label='Carriage type',n.bed=para$n.bed, n.day=para$n.day, timestep=1)
    mosaicdata.car.short$`Carriage type`=factor(mosaicdata.car.short$`Carriage type`, levels = c('sR', 'sr','Sr', 'ss','S'))
    colo_table_filled_long = dat$colo_table_filled_long[[1]]
    mosaicdata.car.long = dataformosaic(data= colo_table_filled_long,label='Carriage type',n.bed=para$n.bed, n.day=para$n.day, timestep=1)
    mosaicdata.car.long$`Carriage type`=factor(mosaicdata.car.long$`Carriage type`, levels = c('sR', 'sr','Sr', 'ss','S'))
    
    sunflower=c('#ee7a12', "#F2A359", "#f8caa0",  "#d6e1d9", "#AAC0AF")
    
    p.car.short=ggplot(mosaicdata.car.short, aes(`Time (days)`, `Bed number`)) + 
        geom_tile(aes(fill = `Carriage type`)) + 
        scale_fill_manual(values=sunflower, name='Carriage type')+
        theme_bw(base_size=base_size) + 
        labs(x = "Observation days", y="Unique patients")+ 
        scale_x_continuous(expand = c(0, 0), breaks= seq(10, para$n.day, 10), labels = as.character(seq(10, para$n.day, 10))) +
        scale_y_continuous(expand = c(0, 0), breaks= 1:10) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.title = element_text(size = base_size+2), 
              legend.text  = element_text(size = base_size+2),
              axis.text = element_text(size = base_size+2, colour = "grey50"), 
              axis.title = element_text(size=base_size+2))+
        geom_segment(data=day1.short, aes(x,y,xend=x.end, yend=y.end), size=0.25, inherit.aes=F, colour='blue')
    
    p.car.long=ggplot(mosaicdata.car.long, aes(`Time (days)`, `Bed number`)) + 
        geom_tile(aes(fill = `Carriage type`)) + 
        scale_fill_manual(values=sunflower,name='Carriage type')+
        theme_bw(base_size=base_size) + 
        labs(x = "Observation days", y="")+ 
        scale_x_continuous(expand = c(0, 0), breaks= seq(10, para$n.day, 10), labels = as.character(seq(10, para$n.day, 10))) +
        scale_y_continuous(expand = c(0, 0), breaks= 1:10) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.title = element_text(size = base_size+2), 
              legend.text  = element_text(size = base_size+2),
              axis.text = element_text(size = base_size+2, colour = "grey50"), 
              axis.title = element_text(size=base_size+2))+
        geom_segment(data=day1.long, aes(x,y,xend=x.end, yend=y.end), size=0.25, inherit.aes=F, colour='blue')
    
    car.mosaic = ggarrange(p.car.short, p.car.long, ncol=2, common.legend = T, legend = 'bottom')
    
    ### PLOT 3 ###
    ##total R per day 
    Rperday = list()
    Rperdayline = list()
    for (t in c('colo_table_filled_short', 'colo_table_filled_long')){
        colo_table_filled= dat[[t]]
        totalRperday = lapply(colo_table_filled, function(x) {
            apply(x, 1, function(x) length(which(x=='sR')))
        })
        Rperdaydata = as.data.frame(totalRperday)
        colnames(Rperdaydata) = 1:ncol(Rperdaydata)
        Rperday[[t]] = Rperdaydata
        Rperdaydata$time = 1:para$n.day
        Rperdaydata.melt = reshape2::melt(Rperdaydata, id.var = c('time'))
        Rperdaydata.melt$dark = 'n'
        Rperdaydata.melt$dark[which(Rperdaydata.melt$variable == 1)] = 'y'
        
        Rperdayline[[t]] = ggplot(Rperdaydata.melt, aes(x = time, y = value, group = variable, alpha = dark, size = dark))+
            geom_line() + 
            labs(x = "Observation days", y = "Number of\nresistance carriers") +
            scale_y_continuous(limit = c(0,10), breaks= 1:10)+
            scale_x_continuous(limit = c(0,para$n.day))+
            scale_alpha_manual(values = c(0.2, 1)) +
            scale_size_manual(values = c(0.2, 1)) +
            theme_bw() +
            guides(size = 'none',
                   alpha = 'none') +
            theme(axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size=base_size+2))
        
        if (t == 'colo_table_filled_long') {Rperdayline[[t]] = Rperdayline[[t]] + theme(axis.title.y = element_blank())}
    }
    
    totalRline = ggarrange(plotlist = Rperdayline, ncol = 2, common.legend = T, legend = 'bottom')
    
    ### PLOT 4 ###
    ##Cumulative R 
    sumdata = list()
    for (t in c('colo_table_filled_short', 'colo_table_filled_long')){
        cs = cumsum(Rperday[[t]])
        cs$time = 1:para$n.day
        cs$dur = t
        sumdata[[t]] = reshape2::melt(cs, id.var = c('time', 'dur'))
    }

    cumsumdata = do.call('rbind.data.frame', sumdata)
    cumsumdata$dark = 'n'
    cumsumdata$dark[which(cumsumdata$variable == 1)] = 'y'
    cumsumdata$group = paste(cumsumdata$variable, cumsumdata$dur)
    sumplot = ggplot(cumsumdata, aes(x = time, y = value, group = group, colour = dur, alpha = dark, size = dark))+
        geom_line() +
        labs(x = "Observation days", y="Cumulative sum of\nresistance carriers")+
        scale_color_manual(values =  c('#D65780', '#93B5C6'), labels = c('Long', 'Short'), name = 'Treatment duration') +
        scale_alpha_manual(values = c(0.1, 1)) +
        scale_size_manual(values = c(0.6, 1)) +
        theme_bw()+
        guides(size = 'none',
               alpha = 'none') + 
        theme(axis.text = element_text(size = base_size+2, colour = "grey50"), 
              axis.title = element_text(size = base_size+2), 
              legend.position = "bottom", 
              legend.background = element_rect(fill=lgdcol,size=0.05),
              legend.text  = element_text(size = base_size+2))
    
    allplots = ggarrange(ggarrange(abx.mosaic, car.mosaic, totalRline, nrow=3), 
                        sumplot, nrow=2, heights = c(3, 1))
    
    return(allplots)
}


