##########################################################################
#######Effect of antibiotic duration on hospitalised patients ############
##################### Show transmission changes    #######################
##########################################################################

# clean working environment
rm(list=ls()) 

# load libraries 
library(ggplot2)
library(ggpubr)
library(reshape)

# load output from model 
load('runs/cocarriage5state_transheat_abs_2021-10-29.Rdata')
d.wide = data.frame(t = unlist(dat['pi_ssr']),
                    i = unlist(dat['prop_R']),
                    d = unlist(dat['meanDur']), 
                    S = unlist(dat['prop S per day']),
                    s = unlist(dat['prop s per day']),
                    sr = unlist(dat['prop sr per day']),
                    Sr = unlist(dat['prop Sr per day']),
                    sR = unlist(dat['prop sR per day']), 
                    p =  unlist(dat['p']))
d = reshape2::melt(d.wide, id.var = c('t', 'i', 'd', 'p'))
d = reshape2::melt(d, id.var = c('d', 'p', 'variable', 'value'))
colnames(d) = c('x', 'p', 'state', 'fill', 'para.vary', 'y' )
d = rbind.data.frame(d[which(d$p == 'pi_ssr' & d$para.vary == 't'),], 
                     d[which(d$p == 'prop_R' & d$para.vary == 'i'),])

counter = 1
plots = list()
for (p in unique(d$p)){
  for (s in c('S', 'Sr', 's', 'sr', 'sR')){
    
    temp.d = d[which(d$p == p & d$state == s),]
    
    plot = ggplot(temp.d, aes(x = x, y = y, fill = fill)) + 
      geom_raster(interpolate=TRUE) + 
      theme_minimal() + 
      scale_fill_gradientn(colours = c('white',  '#f5dea6', '#e39910'),
                           #name = 'Proportion of ward patients with given carriage status',
                           name = '',
                           breaks=c(min(temp.d$fill), median(temp.d$fill), max(temp.d$fill)),
                           labels=c(paste0(round(min(temp.d$fill), 2)*100, '%'),
                                    '',
                                    paste0(round(max(temp.d$fill), 2)*100, '%'))) + 
      theme(text = element_text(size = 12),
            panel.grid.major = element_blank(),
            legend.key.height = unit(0.5, "cm"),
            legend.key.width = unit(0.2, "cm"),
            panel.border = element_blank(),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(), 
            legend.text=element_text(size=7),
            legend.margin=margin(0,0,0,0),
            legend.box.margin=margin(-10,-5,-10,-12)) + 
      guides(fill = guide_colorbar(draw.ulim = FALSE, draw.llim = FALSE, 
                                     ticks = FALSE))
    
    if (counter == 1) {
      plots[[counter]] = plot + 
        labs (y = 'Transmission rate\n', x = '', title = s) +
        theme(legend.position = 'right',
              axis.text.x = element_blank(),
              axis.title.y = element_text(face = 'bold'),
              plot.title = element_text(hjust = 0.5, face= 'bold')) 
    }
    
    if (counter %in% 2:5){
      plots[[counter]] = plot + 
        labs (y = 'Transmission rate', x = '', title = s) +
        theme(legend.position = 'right',
              axis.text = element_blank(), 
              axis.title.y = element_blank(),
              plot.title = element_text(hjust = 0.5, face= 'bold')) 
    }
    
    if (counter == 6){
      plots[[counter]] = plot + 
        labs (y = 'Prevalence of resistance\ncarriers admitted to the ward', x = '') +
        theme(legend.position = 'right', 
              axis.title.y = element_text(face = 'bold')) 
    }
    
    if (counter %in% c(7, 9, 10)){
      plots[[counter]] = plot + 
        labs (y = '', x = '') +
        theme(legend.position = 'right',
              axis.text.y = element_blank(), 
              axis.title.y = element_blank()) 
    }
    
    if (counter == 8 ){
      plots[[counter]] = plot + 
        labs (y = '', x = 'Antibiotic treatment duration (days)') +
        theme(legend.position = 'right',
              axis.text.y = element_blank(), 
              axis.title.y = element_blank()) 
    }
    
    counter = counter + 1
    
  }
}

transheatstates = ggarrange(plotlist=plots, nrow = 2, ncol = 5, widths = c(1.2, 1, 1, 1, 1) )
transheatstates = annotate_figure(transheatstates, top = text_grob('A. Proportion of ward patients with given carriage status', 
                                                     face = "bold", size = 15, hjust = 0.88))
#transheatstates
ggsave('~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/final_main/transheatstate.png',  
       dpi = 500, 
       width = 14, height = 5)


transheat = readRDS('~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/final_main/transheat.rds')
ggarrange(transheatstates, transheat, ncol = 1, nrow = 2, heights = c(1, 1.3))
ggsave('~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/final_main/transheatcomb.png',  
       dpi = 500, 
       width = 12, height = 12)

