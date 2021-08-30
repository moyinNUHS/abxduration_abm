#########################################################################
#######Effect of antibiotic duration on hospitalised patients############
#####################Show one run per scenario ##########################
#########################################################################

dataformosaic <- function(data, label, n.bed = para$n.bed, n.day = para$n.day, timestep=timestep){
    
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

plot_scenario <- function(model, scenario) {
    
    if(model == "simple3state"){
        
        if (scenario == 'esbl') {
            para =  c(common.para$esbl, simple3state.para$esbl)
        } else {
            para =  c(common.para$cpe, simple3state.para$cpe)
        }
        
        #Run model 
        timestep = 1
        n.day = 365
        sdDur = 1
        
        matrixes.short = los_abx_table_varydur(n.bed=para$n.bed, n.day=para$n.day, max.los=para$max.los, 
                                               p.infect=para$p.infect, p.r.day1=para$p.r.day1, cum.r.1=para$cum.r.1, 
                                               meanDur= para$short_dur, timestep=timestep)
        patient.matrix.short = matrixes.short[[1]]
        day1.short = admitdays(patient.matrix.short)
        abx.matrix.short = matrixes.short[[2]]
        los.array.short = summary_los(patient.matrix=patient.matrix.short)
        colo.matrix.short = colo.table(patient.matrix=patient.matrix.short, los=los.array.short, 
                                       prop_R=para$prop_R, prop_S=para$prop_S)
        
        colo_table_filled_short = nextDay(patient.matrix = patient.matrix.short, los.array = los.array.short, 
                                          abx.matrix = abx.matrix.short, colo.matrix = colo.matrix.short, 
                                          bif = para$bif, pi_ssr = para$pi_ssr, repop.s = para$repop.s, mu =para$mu, abx.s = para$abx.s, abx.r = para$abx.r, 
                                          timestep = timestep)
        
        matrixes.long = los_abx_table_varydur(n.bed=para$n.bed, n.day=para$n.day, max.los=para$max.los, 
                                              p.infect=para$p.infect, p.r.day1=para$p.r.day1, cum.r.1=para$cum.r.1, 
                                              meanDur= para$long_dur, timestep=timestep)
        patient.matrix.long = matrixes.long[[1]]
        day1.long = admitdays(patient.matrix.long)
        abx.matrix.long = matrixes.long[[2]]
        los.array.long = summary_los(patient.matrix = patient.matrix.long)
        colo.matrix.long = colo.table(patient.matrix = patient.matrix.long, los = los.array.long, 
                                      prop_R = para$prop_R,prop_S = para$prop_S)
        
        colo_table_filled_long = nextDay(patient.matrix = patient.matrix.long, los.array = los.array.long, 
                                         abx.matrix = abx.matrix.long, colo.matrix = colo.matrix.long, 
                                         bif = para$bif, pi_ssr = para$pi_ssr, repop.s = para$repop.s, mu = para$mu, 
                                         abx.s = para$abx.s, abx.r = para$abx.r,timestep = timestep)
        
        #Plots 
        ####Abx use plots 
        mosaicdata.abx.short = dataformosaic(data=abx.matrix.short ,label='Antibiotic type',n.bed=para$n.bed, n.day=para$n.day, timestep=timestep)
        mosaicdata.abx.short$`Antibiotic type`=as.factor(mosaicdata.abx.short$`Antibiotic type`)
        mosaicdata.abx.long = dataformosaic(data=abx.matrix.long ,label='Antibiotic type',n.bed=para$n.bed, n.day=para$n.day, timestep=timestep)
        mosaicdata.abx.long$`Antibiotic type`=as.factor(mosaicdata.abx.long$`Antibiotic type`)
        
        winteralmond = c('white',"#87C2BE","#5E8E7B")
        lgdcol = rgb(0.185, 0.188, 0.154, alpha = .05)
        base_size = 9
        
        p.abx.short = ggplot(mosaicdata.abx.short, aes(`Time (days)`, `Bed number`)) + 
            geom_tile(aes(fill = `Antibiotic type`)) + 
            scale_fill_manual(values=winteralmond, labels = c("none", "narrow-spectrum", "broad-spectrum"),limits = c("0", "1", "2"), name='Antibiotic type')+
            theme_bw(base_size=base_size) + 
            labs(title = 'Short antibiotic treatment duration',
                 x = "Time (days)", y="Bed number")+ 
            scale_x_continuous(expand = c(0, 0), breaks= c(10, 20, 30)*timestep, labels = c('10','20','30')) +
            scale_y_continuous(expand = c(0, 0), breaks= c(5, 10), labels = c('5','10')) + 
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
                  legend.position = "bottom", 
                  legend.background = element_rect(fill=lgdcol,size=0.05),
                  legend.title = element_text(size = base_size+2), 
                  legend.text  = element_text(size = base_size+2),
                  axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size= base_size+2))+
            geom_segment(data = day1.short, aes(x, y , xend = x.end, yend = y.end), size = 0.15, inherit.aes=F, colour='blue')
        
        p.abx.long = ggplot(mosaicdata.abx.long, aes(`Time (days)`, `Bed number`)) + 
            geom_tile(aes(fill = `Antibiotic type`)) + 
            scale_fill_manual(values=winteralmond, labels = c("none", "narrow-spectrum", "broad-spectrum"),limits = c("0", "1", "2"), name='Antibiotic type')+
            theme_bw(base_size = base_size) + 
            labs(title = 'Long antibiotic treatment duration',
                 x = "Time (days)", y="")+ 
            scale_x_continuous(expand = c(0, 0), breaks= c(10, 20, 30)*timestep, labels = c('10','20','30')) +
            scale_y_continuous(expand = c(0, 0), breaks= c(5, 10), labels = c('5','10')) + 
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
                  legend.position = "bottom", 
                  legend.background = element_rect(fill=lgdcol,size=0.05),
                  legend.title = element_text(size = base_size+2), 
                  legend.text  = element_text(size = base_size+2),
                  axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size=base_size+2))+
            geom_segment(data=day1.long, aes(x,y,xend=x.end, yend=y.end), size=0.25, inherit.aes=F, colour='blue')
        
        abx.mosaic=ggarrange(p.abx.short, p.abx.long, ncol=2, common.legend = T, legend = 'bottom')
        
        ##Carriage mosaic 
        mosaicdata.car.short = dataformosaic(data=colo_table_filled_short ,label='Carriage type',n.bed=para$n.bed, n.day=para$n.day, timestep=timestep)
        mosaicdata.car.long = dataformosaic(data= colo_table_filled_long ,label='Carriage type',n.bed=para$n.bed, n.day=para$n.day, timestep=timestep)
        
        sunflower = c('#F2A359',"#AAC0AF","#d6e1d9")
        
        p.car.short = ggplot(mosaicdata.car.short, aes(`Time (days)`, `Bed number`)) + 
            geom_tile(aes(fill = `Carriage type`)) + 
            scale_fill_manual(values=sunflower,labels=c('R','S','ss'),limits = c('R','S','ss'),name='Carriage type')+
            theme_bw(base_size=base_size) + 
            labs(x = "Time (days)", y="Bed number")+ 
            scale_x_continuous(expand = c(0, 0), breaks= c(10, 20, 30)*timestep, labels = c('10','20','30')) +
            scale_y_continuous(expand = c(0, 0), breaks= c(5, 10), labels = c('5','10')) +
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
            scale_fill_manual(values=sunflower,labels=c('R','S','ss'),limits = c('R','S','ss'),name='Carriage type')+
            theme_bw(base_size=base_size) + 
            labs(x = "Time (days)", y="Bed number")+ 
            scale_x_continuous(expand = c(0, 0), breaks= c(10, 20, 30)*timestep, labels = c('10','20','30')) +
            scale_y_continuous(expand = c(0, 0), breaks= c(5, 10), labels = c('5','10')) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
                  legend.position = "bottom", 
                  legend.background = element_rect(fill=lgdcol,size=0.05),
                  legend.title = element_text(size = base_size+2), 
                  legend.text  = element_text(size = base_size+2),
                  axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size=base_size+2))+
            geom_segment(data=day1.long, aes(x,y,xend=x.end, yend=y.end), size=0.25, inherit.aes=F, colour='blue')
        
        car.mosaic=ggarrange(p.car.short, p.car.long, ncol=2, common.legend = T, legend = 'bottom')
        
        ##total R per day 
        totalRperday.short=apply(colo_table_filled_short, 1, function(x) length(which(x=='R')))
        totalRperday.short.avg= rowMeans(matrix(totalRperday.short, ncol=timestep, byrow=T))
        Rperdaydata.short=data.frame(`Time (days)`= 1:para$n.day, 
                                     `Total R per day`= totalRperday.short.avg)
        Rperdayline.short=ggplot(Rperdaydata.short, aes(x=Time..days., y=Total.R.per.day))+
            geom_point(colour= 'grey50', size=0.1)+ 
            geom_line(colour= 'grey50') + 
            labs(x = "Time (days)", y="Number of R carrier\nper day")+
            scale_y_continuous(limit = c(0,10), breaks= c(5, 10))+
            theme_bw()+
            theme(axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size=base_size+2))
        
        totalRperday.long=apply(colo_table_filled_long, 1, function(x) length(which(x=='R')))
        totalRperday.long.avg= rowMeans(matrix(totalRperday.long, ncol=timestep, byrow=T))
        Rperdaydata.long=data.frame(`Time (days)`= 1:para$n.day, 
                                    `Total R per day`= totalRperday.long.avg)
        Rperdayline.long=ggplot(Rperdaydata.long, aes(x=Time..days., y=Total.R.per.day))+
            geom_point(colour= 'grey50', size=0.1)+ 
            geom_line(colour= 'grey50') + 
            labs(x = "Time (days)", y="")+
            scale_y_continuous(limit = c(0,10), breaks= c(5, 10))+
            theme_bw()+
            theme(axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size=base_size+2))
        
        totalRline=ggarrange(Rperdayline.short, Rperdayline.long, ncol=2, common.legend = T, legend = 'bottom')
        
        ##Cumulative R 
        cumsum.short=cumsum(totalRperday.short.avg)
        cumsumdata.short=data.frame(`Time (days)`= 1:para$n.day, 
                                    cumsum= cumsum.short)
        cumsumdata.short$dur=rep('Short',nrow(cumsumdata.short))
        cumsum.long=cumsum(totalRperday.long.avg)
        cumsumdata.long=data.frame(`Time (days)`= 1:para$n.day, 
                                   cumsum= cumsum.long)
        cumsumdata.long$dur=rep('Long',nrow(cumsumdata.long))
        cumsumdata=rbind.data.frame(cumsumdata.short,cumsumdata.long)
        sumplot=ggplot(cumsumdata, aes(x=Time..days., y=cumsum, colour = dur))+
            geom_line()+
            labs(x = "Time (days)", y="Cumulative sum of R \ncarriers")+
            scale_color_manual(values =  c('#A0495B', '#0096BC'), name= 'Treatment duration')+
            theme_bw()+
            theme(axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size=base_size+2), 
                  legend.position = "bottom", 
                  legend.background = element_rect(fill=lgdcol,size=0.05),
                  legend.text  = element_text(size = base_size+2))
        
        (allplots= ggarrange(ggarrange(abx.mosaic, car.mosaic, totalRline, nrow=3), 
                             sumplot, nrow=2, heights = c(3, 1)))
        
    }else if(model == "cocarriage5state"){
        
        if (scenario == 'esbl') {
            para =  c(common.para$esbl, cocarriage5state.para$esbl)
        } else {
            para = c(common.para$cpe, cocarriage5state.para$cpe)
        }
        
        #Run model 
        timestep = 1
        n.day = 365
        sdDur = 1
        
        matrixes = los_abx_table_varydur(n.bed = para$n.bed, n.day = para$n.day, max.los = para$max.los, 
                                         p.infect = para$p.infect, p.r.day1 = para$p.r.day1, p.r.after = para$p.r.after, cum.r.1 = para$cum.r.1, 
                                         meanDur= para$short_dur, timestep = timestep)
        patient.matrix.short = matrixes[[1]]
        day1.short = admitdays(patient.matrix.short)
        
        abx.matrix.short = matrixes[[2]]
        los.array.short = summary_los(patient.matrix=patient.matrix.short)
        colo.matrix.short = colo.table(patient.matrix=patient.matrix.short, los=los.array.short, 
                                       prop_R=para$prop_R, prop_r=para$prop_r, prop_Sr=para$prop_Sr, prop_S=para$prop_S)
        
        colo_table_filled_short = nextDay(patient.matrix=patient.matrix.short, abx.matrix=abx.matrix.short, colo.matrix=colo.matrix.short, 
                                          pi_ssr=para$pi_ssr, bif=para$bif, mu=para$mu, repop.r=para$repop.r,
                                          repop.s=para$repop.s, abx.r=para$abx.r, abx.s=para$abx.s, timestep=timestep)
        
        matrixes = los_abx_table_varydur(n.bed=para$n.bed, n.day=para$n.day, max.los=para$max.los, 
                                         p.infect=para$p.infect, p.r.day1=para$p.r.day1, cum.r.1=para$cum.r.1, 
                                         meanDur= para$long_dur, timestep=timestep)
        patient.matrix.long=matrixes[[1]]
        day1.long = admitdays(patient.matrix.long)
        
        abx.matrix.long=matrixes[[2]]
        los.array.long = summary_los(patient.matrix=patient.matrix.long)
        colo.matrix.long = colo.table(patient.matrix=patient.matrix.long, los=los.array.long, 
                                      prop_R=para$prop_R, prop_r=para$prop_r, prop_Sr=para$prop_Sr, prop_S=para$prop_S)
        
        colo_table_filled_long = nextDay(patient.matrix=patient.matrix.long, abx.matrix=abx.matrix.long, colo.matrix=colo.matrix.long, 
                                         pi_ssr=para$pi_ssr, bif=para$bif, mu=para$mu, repop.r=para$repop.r,
                                         repop.s=para$repop.s, abx.r=para$abx.r, abx.s=para$abx.s, timestep=timestep)
        
        #Plots 
        ####Abx use plots 
        mosaicdata.abx.short = dataformosaic(data=abx.matrix.short ,label='Antibiotic type',n.bed=para$n.bed, n.day=para$n.day, timestep=timestep)
        mosaicdata.abx.short$`Antibiotic type`=as.factor(mosaicdata.abx.short$`Antibiotic type`)
        mosaicdata.abx.long = dataformosaic(data=abx.matrix.long ,label='Antibiotic type',n.bed=para$n.bed, n.day=para$n.day, timestep=timestep)
        mosaicdata.abx.long$`Antibiotic type`=as.factor(mosaicdata.abx.long$`Antibiotic type`)
        
        winteralmond = c('white',"#87C2BE","#5E8E7B")
        lgdcol = rgb(0.185, 0.188, 0.154, alpha = .05)
        base_size = 9
        
        p.abx.short = ggplot(mosaicdata.abx.short, aes(`Time (days)`, `Bed number`)) + 
            geom_tile(aes(fill = `Antibiotic type`)) + 
            scale_fill_manual(values=winteralmond,labels = c("None", "Third-generation cephalosporin", "Carbapenem"), name='Antibiotic type')+
            theme_bw(base_size=base_size) + 
            labs(title= 'Short antibiotic treatment duration',
                 x = "Observation days", y = "Bed number")+ 
            scale_x_continuous(expand = c(0, 0), breaks= c(10,20,30)*timestep, labels = c('10','20','30')) +
            scale_y_continuous(expand = c(0, 0), breaks= c(5, 10), labels = c('5','10')) +
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
            scale_fill_manual(values=winteralmond,labels = c("None", "Third-generation cephalosporin", "Carbapenem"), name='Antibiotic type')+
            theme_bw(base_size=base_size) + 
            labs(title= 'Long antibiotic treatment duration',
                 x = "Observation days", y="")+ 
            scale_x_continuous(expand = c(0, 0), breaks= c(10,20,30)*timestep, labels = c('10','20','30')) +
            scale_y_continuous(expand = c(0, 0), breaks= c(5, 10), labels = c('5','10')) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
                  legend.position = "bottom", 
                  legend.background = element_rect(fill=lgdcol,size=0.05),
                  legend.title = element_text(size = base_size+2), 
                  legend.text  = element_text(size = base_size+2),
                  axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size=base_size+2))+
            geom_segment(data=day1.long, aes(x,y,xend=x.end, yend=y.end), size=0.25, inherit.aes=F, colour='blue')
        
        
        abx.mosaic=ggarrange(p.abx.short, p.abx.long, ncol=2, common.legend = T, legend = 'bottom')
        
        ##Carriage mosaic 
        mosaicdata.car.short = dataformosaic(data=colo_table_filled_short,label='Carriage type',n.bed=para$n.bed, n.day=para$n.day, timestep=timestep)
        mosaicdata.car.short$`Carriage type`=factor(mosaicdata.car.short$`Carriage type`, levels = c('sR', 'sr','Sr', 'ss','S'))
        mosaicdata.car.long = dataformosaic(data= colo_table_filled_long,label='Carriage type',n.bed=para$n.bed, n.day=para$n.day, timestep=timestep)
        mosaicdata.car.long$`Carriage type`=factor(mosaicdata.car.long$`Carriage type`, levels = c('sR', 'sr','Sr', 'ss','S'))
        
        sunflower=c('#ee7a12', "#F2A359", "#f8caa0",  "#d6e1d9", "#AAC0AF")
        
        p.car.short=ggplot(mosaicdata.car.short, aes(`Time (days)`, `Bed number`)) + 
            geom_tile(aes(fill = `Carriage type`)) + 
            scale_fill_manual(values=sunflower, name='Carriage type')+
            theme_bw(base_size=base_size) + 
            labs(x = "Observation days", y="Bed number")+ 
            scale_x_continuous(expand = c(0, 0), breaks= c(10,20,30)*timestep, labels = c('10','20','30')) +
            scale_y_continuous(expand = c(0, 0), breaks= c(5, 10), labels = c('5','10')) +
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
            scale_x_continuous(expand = c(0, 0), breaks= c(10, 20, 30)*timestep, labels = c('10','20','30')) +
            scale_y_continuous(expand = c(0, 0), breaks= c(5, 10), labels = c('5','10')) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
                  legend.position = "bottom", 
                  legend.background = element_rect(fill=lgdcol,size=0.05),
                  legend.title = element_text(size = base_size+2), 
                  legend.text  = element_text(size = base_size+2),
                  axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size=base_size+2))+
            geom_segment(data=day1.long, aes(x,y,xend=x.end, yend=y.end), size=0.25, inherit.aes=F, colour='blue')
        
        car.mosaic=ggarrange(p.car.short, p.car.long, ncol=2, common.legend = T, legend = 'bottom')
        
        ##total R per day 
        totalRperday.short=apply(colo_table_filled_short, 1, function(x) length(which(x=='sR')))
        totalRperday.short.avg= rowMeans(matrix(totalRperday.short, ncol=timestep, byrow=T))
        Rperdaydata.short=data.frame(`Time (days)`= 1:para$n.day, 
                                     `Total R per day`= totalRperday.short.avg)
        Rperdayline.short=ggplot(Rperdaydata.short, aes(x=Time..days., y=Total.R.per.day))+
            geom_point(colour= 'grey50', size=0.1)+ 
            geom_line(colour= 'grey50') + 
            labs(x = "Observation days", y="Number of\nresistance carriers")+
            scale_y_continuous(limit = c(0,10), breaks= c(5, 10))+
            theme_bw()+
            theme(axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size=base_size+2))
        
        totalRperday.long=apply(colo_table_filled_long, 1, function(x) length(which(x=='sR')))
        totalRperday.long.avg= rowMeans(matrix(totalRperday.long, ncol=timestep, byrow=T))
        Rperdaydata.long=data.frame(`Time (days)`= 1:para$n.day, 
                                    `Total R per day`= totalRperday.long.avg)
        Rperdayline.long=ggplot(Rperdaydata.long, aes(x=Time..days., y=Total.R.per.day))+
            geom_point(colour= 'grey50', size=0.1)+ 
            geom_line(colour= 'grey50') + 
            labs(x = "Observation days", y="")+
            scale_y_continuous(limit = c(0,10), breaks= c(5, 10))+
            theme_bw()+
            theme(axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size=base_size+2))
        
        totalRline=ggarrange(Rperdayline.short, Rperdayline.long, ncol=2, common.legend = T, legend = 'bottom')
        
        ##Cumulative R 
        cumsum.short=cumsum(totalRperday.short.avg)
        cumsumdata.short=data.frame(`Time (days)`= 1:para$n.day, 
                                    cumsum = cumsum.short)
        cumsumdata.short$dur=rep('Short',nrow(cumsumdata.short))
        cumsum.long=cumsum(totalRperday.long.avg)
        cumsumdata.long=data.frame(`Time (days)`= 1:para$n.day, 
                                   cumsum= cumsum.long)
        cumsumdata.long$dur=rep('Long',nrow(cumsumdata.long))
        cumsumdata=rbind.data.frame(cumsumdata.short,cumsumdata.long)
        sumplot=ggplot(cumsumdata, aes(x=Time..days., y=cumsum, colour = dur))+
            geom_line()+
            labs(x = "Observation days", y="Cumulative sum of\nresistance carriers")+
            scale_color_manual(values =  c('#D65780', '#93B5C6'), name='Treatment duration')+
            theme_bw()+
            theme(axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size=base_size+2), 
                  legend.position = "bottom", 
                  legend.background = element_rect(fill=lgdcol,size=0.05),
                  legend.text  = element_text(size = base_size+2))
        
        (allplots= ggarrange(ggarrange(abx.mosaic, car.mosaic, totalRline, nrow=3), 
                             sumplot, nrow=2, heights = c(3, 1)))
        
        
    }else if(model == "populationgrowth"){
        
        if (scenario == 'esbl') {
            para = c(common.para$esbl, populationgrowth.para$esbl)
        } else {
            para = c(common.para$cpe, populationgrowth.para$cpe)
        }
        
        timestep=1
        
        matrixes = los_abx_table_varydur(n.bed=para$n.bed, n.day=para$n.day, max.los=para$max.los, 
                                         p.infect=para$p.infect, p.r.day1=para$p.r.day1, cum.r.1=para$cum.r.1, 
                                         meanDur= para$short_dur, timestep=timestep)
        patient.matrix.short=matrixes[[1]]
        day1.short= admitdays(patient.matrix.short)
        abx.matrix.short=matrixes[[2]]
        los.array.short = summary_los(patient.matrix=patient.matrix.short)
        colo.matrix.short = colo.table(patient.matrix=patient.matrix.short, los.array=los.array.short, total_prop=para$total_prop, prop_R=para$prop_R,r_thres=para$r_thres, K=para$K)
        colo_table_filled_short = nextDay(patient.matrix=patient.matrix.short, los.array=los.array.short, abx.matrix=abx.matrix.short, colo.matrix=colo.matrix.short, 
                                          pi_ssr=para$pi_ssr, total_prop = para$total_prop, K=para$K, r_growth=para$r_growth, r_thres=para$r_thres, r_trans=para$r_trans,s_growth=para$s_growth,
                                          abx.s=para$abx.s, abx.r=para$abx.r, timestep=timestep)[[2]]
        
        matrixes = los_abx_table_varydur(n.bed = para$n.bed, n.day = para$n.day, max.los = para$max.los, 
                                         p.infect = para$p.infect, p.r.day1 = para$p.r.day1, cum.r.1 = para$cum.r.1, 
                                         meanDur = para$long_dur, timestep = timestep)
        patient.matrix.long=matrixes[[1]]
        day1.long= admitdays(patient.matrix.long)
        abx.matrix.long=matrixes[[2]]
        los.array.long = summary_los(patient.matrix=patient.matrix.long)
        colo.matrix.long = colo.table(patient.matrix=patient.matrix.long, los.array=los.array.long, total_prop=para$total_prop, prop_R=para$prop_R,r_thres=para$r_thres, K=para$K)
        colo_table_filled_long = nextDay(patient.matrix=patient.matrix.long, los.array=los.array.long, abx.matrix=abx.matrix.long, colo.matrix=colo.matrix.long, 
                                         pi_ssr=para$pi_ssr, total_prop = para$total_prop, K=para$K, r_growth=para$r_growth, r_thres=para$r_thres, r_trans=para$r_trans,s_growth=para$s_growth,
                                         abx.s=para$abx.s, abx.r=para$abx.r, timestep=timestep)[[2]]
        
        #Plots 
        ####Abx use plots 
        mosaicdata.abx.short = dataformosaic(data=abx.matrix.short ,label='Antibiotic type',n.bed=para$n.bed, n.day=para$n.day, timestep=timestep)
        mosaicdata.abx.short$`Antibiotic type`=as.factor(mosaicdata.abx.short$`Antibiotic type`)
        mosaicdata.abx.long = dataformosaic(data=abx.matrix.long ,label='Antibiotic type',n.bed=para$n.bed, n.day=para$n.day, timestep=timestep)
        mosaicdata.abx.long$`Antibiotic type`=as.factor(mosaicdata.abx.long$`Antibiotic type`)
        
        winteralmond = c('white',"#87C2BE","#5E8E7B")
        lgdcol = rgb(0.185, 0.188, 0.154, alpha = .05)
        base_size = 9
        
        p.abx.short=ggplot(mosaicdata.abx.short, aes(`Time (days)`, `Bed number`)) + 
            geom_tile(aes(fill = `Antibiotic type`)) + 
            scale_fill_manual(values=winteralmond,labels = c("none", "narrow-spectrum", "broad-spectrum"), name='Antibiotic type')+
            theme_bw(base_size=base_size) + 
            labs(title= 'Short antibiotic treatment duration',
                 x = "Time (days)", y="Bed number")+ 
            scale_x_continuous(expand = c(0, 0), breaks= c(10, 20, 30)) +
            scale_y_continuous(expand = c(0, 0), breaks= c(5, 10), labels = c('5','10')) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
                  legend.position = "bottom", 
                  legend.background = element_rect(fill=lgdcol,size=0.05),
                  legend.title = element_text(size = base_size+2), 
                  legend.text  = element_text(size = base_size+2),
                  axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size=base_size+2))+
            geom_segment(data=day1.short, aes(x,y,xend=x.end, yend=y.end), size=0.25, inherit.aes=F, colour='blue')
        
        p.abx.long=ggplot(mosaicdata.abx.long, aes(`Time (days)`, `Bed number`)) + 
            geom_tile(aes(fill = `Antibiotic type`)) + 
            scale_fill_manual(values=winteralmond,labels = c("none", "narrow-spectrum", "broad-spectrum"), name='Antibiotic type')+
            theme_bw(base_size=base_size) + 
            labs(title= 'Long antibiotic treatment duration',
                 x = "Time (days)", y="")+ 
            scale_x_continuous(expand = c(0, 0), breaks= c(10, 20, 30)) +
            scale_y_continuous(expand = c(0, 0), breaks= c(5, 10), labels = c('5','10')) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
                  legend.position = "bottom", 
                  legend.background = element_rect(fill=lgdcol,size=0.05),
                  legend.title = element_text(size = base_size+2), 
                  legend.text  = element_text(size = base_size+2),
                  axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size=base_size+2))+
            geom_segment(data=day1.long, aes(x,y,xend=x.end, yend=y.end), size=0.25, inherit.aes=F, colour='blue')
        
        abx.mosaic=ggarrange(p.abx.short, p.abx.long, ncol=2, common.legend = T, legend = 'bottom')
        
        ##Carriage mosaic 
        mosaicdata.car.short = dataformosaic(data=colo_table_filled_short,label='Number of R',n.bed=para$n.bed, n.day=para$n.day, timestep=timestep)
        mosaicdata.car.short$abovethreshold= as.factor(mosaicdata.car.short$`Number of R`>= r_thres)
        
        mosaicdata.car.long = dataformosaic(data=colo_table_filled_long,label='Number of R',n.bed=para$n.bed, n.day=para$n.day, timestep=timestep)
        mosaicdata.car.long$abovethreshold= as.factor(mosaicdata.car.long$`Number of R`>= r_thres)
        
        sunflower=c("#d6e1d9",'#ee7a12')
        
        p.car.short=ggplot(mosaicdata.car.short, aes(`Time (days)`, `Bed number`)) + 
            geom_tile(aes(fill = abovethreshold)) + 
            scale_fill_manual(values=sunflower, name='Carriage type', 
                              labels=c('Below threshold for R transmission','Above threshold for R transmission'))+
            theme_bw(base_size=base_size) + 
            labs(x = "Time (days)", y="Bed number")+ 
            scale_x_continuous(expand = c(0, 0), breaks= c(10, 20, 30)) +
            scale_y_continuous(expand = c(0, 0), breaks= c(5, 10), labels = c('5','10')) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
                  legend.position = "bottom", 
                  legend.background = element_rect(fill=lgdcol,size=0.05),
                  legend.title = element_text(size = base_size+2), 
                  legend.text  = element_text(size = base_size+2),
                  axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size=base_size+2))+
            geom_segment(data=day1.short, aes(x,y,xend=x.end, yend=y.end), size=0.25, inherit.aes=F, colour='blue')
        
        p.car.long=ggplot(mosaicdata.car.long, aes(`Time (days)`, `Bed number`)) + 
            geom_tile(aes(fill = abovethreshold)) + 
            scale_fill_manual(values=sunflower, name='Carriage type', 
                              labels=c('Below threshold for R transmission','Above threshold for R transmission'))+
            theme_bw(base_size=base_size) + 
            labs(x = "Time (days)", y="Bed number")+ 
            scale_x_continuous(expand = c(0, 0), breaks= c(10, 20, 30)) +
            scale_y_continuous(expand = c(0, 0), breaks= c(5, 10), labels = c('5','10')) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size=14),
                  legend.position = "bottom", 
                  legend.background = element_rect(fill=lgdcol,size=0.05),
                  legend.title = element_text(size = base_size+2), 
                  legend.text  = element_text(size = base_size+2),
                  axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size=base_size+2))+
            geom_segment(data=day1.long, aes(x,y,xend=x.end, yend=y.end), size=0.25, inherit.aes=F, colour='blue')
        
        car.mosaic=ggarrange(p.car.short, p.car.long, ncol=2, common.legend = T, legend = 'bottom')
        
        ##total R per day 
        totalRperday.short= rowSums(colo_table_filled_short >= r_thres)
        totalRperday.short.avg= rowMeans(matrix(totalRperday.short, ncol=timestep, byrow=T))
        Rperdaydata.short=data.frame(`Time (days)`= 1:para$n.day, 
                                     `Total R per day`= totalRperday.short.avg)
        Rperdayline.short=ggplot(Rperdaydata.short, aes(x=Time..days., y=Total.R.per.day))+
            geom_point(colour= 'grey50', size=0.1)+ 
            geom_line(colour= 'grey50') + 
            labs(x = "Time (days)", y="Number of R carriers\nper day")+
            scale_y_continuous(limit = c(0,10), breaks= c(5, 10))+
            theme_bw()+
            theme(axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size=base_size+2))
        
        totalRperday.long=rowSums(colo_table_filled_long >= r_thres)
        totalRperday.long.avg= rowMeans(matrix(totalRperday.long, ncol=timestep, byrow=T))
        Rperdaydata.long=data.frame(`Time (days)`= 1:para$n.day, 
                                    `Total R per day`= totalRperday.long.avg)
        Rperdayline.long=ggplot(Rperdaydata.long, aes(x=Time..days., y=Total.R.per.day))+
            geom_point(colour= 'grey50', size=0.1)+ 
            geom_line(colour= 'grey50') + 
            labs(x = "Time (days)", y="")+
            scale_y_continuous(limit = c(0,10), breaks= c(5, 10))+
            theme_bw()+
            theme(axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size=base_size+2))
        
        totalRline=ggarrange(Rperdayline.short, Rperdayline.long, ncol=2, common.legend = T, legend = 'bottom')
        
        ##Cumulative R 
        cumsum.short=cumsum(totalRperday.short.avg)
        cumsumdata.short=data.frame(`Time (days)`= 1:para$n.day, 
                                    cumsum= cumsum.short)
        cumsumdata.short$`Treatment duration`=rep('Short',nrow(cumsumdata.short))
        cumsum.long=cumsum(totalRperday.long.avg)
        cumsumdata.long=data.frame(`Time (days)`= 1:para$n.day, 
                                   cumsum= cumsum.long)
        cumsumdata.long$`Treatment duration`=rep('Long',nrow(cumsumdata.long))
        cumsumdata=rbind.data.frame(cumsumdata.short,cumsumdata.long)
        sumplot=ggplot(cumsumdata, aes(x=Time..days., y=cumsum, colour = `Treatment duration`))+
            geom_line()+
            labs(x = "Time (days)", y="Cumulative sum of R\ncarriers")+
            scale_color_manual(values =  c('#A0495B', '#0096BC'), name='Treatment duration' )+
            theme_bw()+
            theme(axis.text = element_text(size = base_size+2, colour = "grey50"), 
                  axis.title = element_text(size=base_size+2), 
                  legend.position = "bottom", 
                  legend.background = element_rect(fill=lgdcol,size=0.05),
                  legend.text  = element_text(size = base_size+2))
        
        (allplots = ggarrange(ggarrange(abx.mosaic, car.mosaic, totalRline, nrow=3), 
                              sumplot, nrow=2, heights = c(3, 1)))
    }
    return(allplots)
}


