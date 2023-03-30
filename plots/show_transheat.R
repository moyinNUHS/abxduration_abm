##########################################################################
#######Effect of antibiotic duration on hospitalised patients ############
##################### Show transmission changes    #######################
##########################################################################

# clean working environment
rm(list=ls()) 

# load libraries 
library(ggplot2)
library(ggpubr)
library(scales)
library(reshape)

outcomes = c('Treated individuals',
             'Overall ward population', 
             'Non-resistance carriers at admission')

title = 'Effect of longer antibiotic duration on proportion of resistance carriers per day'

highcol = '#CF4D6F'
midcol = '#F1F7ED'
lowcol = '#006989'

eff.lab = 'Antibiotic active\nagainst susceptible\nand resistant organisms'
ineff.lab = 'Antibiotic active\nonly against\nsusceptible organisms'

# run use simple model 
files = c('runs/simple3state_transheat_absdiffnotzero2021-10-30.Rdata', 
          'runs/simple3state_transheat_absdiffzero2021-10-30.Rdata')
d = list()
for (f in files){
  load(f)
  lab = ifelse(length(grep('notzero', f)) == 1, eff.lab, ineff.lab)
  d1 = data.frame(t = unlist(dat['pi_ssr']),
                  i = unlist(dat['prop_R']),
                  f = unlist(dat['totalR_diff']),
                  abx = lab, 
                  outcomes = outcomes[[2]])
  d2 = data.frame(t = unlist(dat['pi_ssr']),
                  i = unlist(dat['prop_R']),
                  f = unlist(dat['newR_diff']),
                  abx = lab, 
                  outcomes = outcomes[[3]])
  d3 = data.frame(t = unlist(dat['pi_ssr']),
                  i = unlist(dat['prop_R']),
                  f = unlist(dat['totalRtreated_diff']),
                  abx = lab, 
                  outcomes = outcomes[[1]])
  d[[f]] = rbind.data.frame(d1, d2, d3)
}

dall = do.call('rbind', d)
dall$outcomes = as.factor(dall$outcomes)
dall$outcomes = factor(dall$outcomes, levels = outcomes)
dall$abx = as.factor(dall$abx)
dall$abx = factor(dall$abx, levels = c(eff.lab, ineff.lab))

lim = max(dall$f)
p = ggplot(dall, aes(x = t, y = i, fill = f)) + 
      labs(x = 'Transmission rate', y = 'Prevalence of resistance carriers on admission') +
      geom_raster(interpolate=TRUE) + 
      theme_minimal() + 
      scale_fill_gradientn(colors = c(lowcol,  midcol, highcol), 
                           values = rescale(c(-lim, 0, lim)),
                           breaks=c(-lim, 0, lim),
                           limit = c(-lim, lim),
                           name = '',
                           labels=c(paste0('More\nresistance\ncarriers\nin the\nshort ward\nthan the\nlong ward\n(-',round(lim,2)*100,'%)'),
                                    'No\ndifference',
                                    paste0('More\nresistance\ncarriers\nin the\nlong ward\nthan the\nshort ward\n(',round(lim,2)*100,'%)'))) +
      facet_grid(abx ~ outcomes, scales="free", switch = "y") +
      theme(text = element_text(size = 13), 
            legend.position = 'right',
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            legend.key.width = unit(0.2,'cm'),
            axis.title = element_text(face = 'bold'),
            legend.key.height = unit(2,'cm'),
            strip.placement = "outside",
            plot.margin = margin(1, 0, 1, 0),
            strip.text = element_text(face = 'bold'),
            panel.background = element_blank()) 

transheat = annotate_figure(p, 
                            top = text_grob('Difference in proportion of resistance carriers between the long and short wards', 
                                            face = "bold", size = 15, hjust = 0.6))

transheat
saveRDS(transheat, file = "~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/final_main/transheat.rds")

ggsave('~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/final_main/transheat.png',  
       dpi = 500, 
       width = 11, height = 7)

