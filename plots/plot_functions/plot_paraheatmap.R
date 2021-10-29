substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# get ranking position of each parameter, add NA if a parameter is not in a particular model
getposition <- function(lhs.filename, labs = labs.df, outcome.type){
  
  d = get(load(lhs.filename))  
  
  if(length(grep('treated', lhs.filename)) == 0) { # get res.names
    source('models/get_output_absdiff_simple3state.R')
    
    if(length(d$prcc) == 6) {
      res.names = res.names[-grep('treated', res.names)]
    }
    
  } else {
    source('models/get_output_treated_simple3state.R')
  }
  
  
  outcome.type.label = which(res.names == outcome.type)
  
  prcc = d$prcc[[outcome.type.label]]$PRCC
  prcc$ranking = NA
  prcc$ranking[which(prcc$original<0)] = scales::rescale(prcc$original[which(prcc$original<0)], 
                                                         to=c(0,0.5)) # scale to 0 to 0.5 if less than 0 
  prcc$ranking[which(prcc$original>0)] = scales::rescale(prcc$original[prcc$original>0], 
                                                         to=c(0.5,1))
  prcc$ranking[which(prcc$`min. c.i.` < 0 & prcc$`max. c.i.`> 0)] = 0.5 #those with CI crossing 0 given 0.5
  
  prcc = prcc[order(prcc$ranking),] #arrange in ranking order
  # low = which(prcc$`min. c.i.` < 0 & prcc$`max. c.i.`< 0)
  # high= which(prcc$`min. c.i.` > 0 & prcc$`max. c.i.`> 0)
  # none= setdiff(1:nrow(prcc),c(low,high))
  # prcc$ranking[low]=prcc$`min. c.i.`[low]
  # prcc$ranking[high]=prcc$`max. c.i.`[high]
  # prcc$ranking[none]=0
  
  ranked.values = prcc$ranking
  names(ranked.values) = rownames(prcc)
  
  missing.para = parameters[!parameters %in% names(ranked.values)]
  to.add.missing.para = rep(NA, length(missing.para))
  names(to.add.missing.para) = missing.para
  
  combine.ranked.values = c(ranked.values, to.add.missing.para)
  
  out = combine.ranked.values[match(levels(parameters), names(combine.ranked.values))]
  
  if (length(grep('treat',outcome.type)) ==1) {outcome.type.lab = 'treated'}
  if (length(grep('totalR_diff',outcome.type)) ==1) {outcome.type.lab = 'totalR'}
  if (length(grep('newR_diff',outcome.type)) ==1) {outcome.type.lab = 'newR'}
  
  names(out) = rev(param.labels)
  name = paste(substr(sub("\\_.*", "", lhs.filename),6,30), ifelse(length(grep('notzero', lhs.filename) == 1), 'notzero', 'zero'), outcome.type.lab)
  
  return(list(out, name))
}

# clean data for plot 
clean_ranks_data <- function(treated, newR, totalR){

  
  # combine the data 
  d = list(treated, totalR, newR)
  
  ## arrange the data
  d.list = lapply(d, function(x){
    
   dat.list = lapply(x, `[[`, 1)
   names(dat.list) = lapply(x, `[[`, 2)
   as.data.frame(do.call('cbind', dat.list))
    
  })
  
  ranks.df.wide = as.data.frame(do.call('cbind', d.list))
  growth.s.gc = ranks.df.wide[rownames(ranks.df.wide) %in% c('c_s', 'g_s'),]
  growth.s = colSums(growth.s.gc, na.rm = T)
  ranks.df.wide[which(rownames(ranks.df.wide) == 'c_s'),] = growth.s
  ranks.df.wide = ranks.df.wide[-which(rownames(ranks.df.wide) == 'g_s'),]
  ranks.df.wide$parameters = rownames(ranks.df.wide)

  ranks.df.wide$empty = NA # create empty first column to place parameter names 
  ranks.df.wide$empty2 = NA
  ranks.df.wide$empty3 = NA
  ranks.df = reshape2::melt(ranks.df.wide, id.var = c('parameters'))
  ranks.df$variable = as.factor(ranks.df$variable)
  models = c('simple3state', 'cocarriage5state', 'populationgrowth')
  ranks.df$variable = factor(ranks.df$variable, 
                             levels = c('empty', 'empty2', 'empty3',
                                        paste0(models, ' zero ', 'treated'), 
                                        paste0(models, ' notzero ', 'treated'),
                                        paste0(models, ' zero ', 'totalR'),
                                        paste0(models, ' notzero ', 'totalR'),
                                        paste0(models, ' zero ', 'newR'), 
                                        paste0(models, ' notzero ', 'newR')))
  ranks.df$value[intersect(which(ranks.df$parameters == 'alpha_r'), grep(' zero', ranks.df$variable))] = NA # in effective abx not available scenario, killing rate of R not relevant
  ranks.df$value[intersect(which(ranks.df$parameters == 'omega_day1.r'), grep(' zero', ranks.df$variable))] = NA # proportion of different abx also not relevant 
  ranks.df$value[intersect(which(ranks.df$parameters == 'omega_after.r'), grep(' zero', ranks.df$variable))] = NA 
  ranks.df$value[intersect(which(ranks.df$parameters == 'omega_day1'), grep(' treated', ranks.df$variable))] = NA # in treated patients, % presribed antibiotics not relevant 
  ranks.df$parameters = as.factor(ranks.df$parameters)
  ranks.df$parameters = factor(ranks.df$parameters, 
                               levels = rev(param.labels[-which(param.labels == 'g_s')]))
  
  #remove parameters that need not be shown 
  to.remove = c('t_short', 't_long', 'p_S', 'p_Sr', 'p_r', 'n', 'l')
  plot.df = ranks.df[-which(ranks.df$parameters %in% to.remove),]
  
  return(plot.df)
}


# plot function 
plot_paraheatmap <- function (plot.df) { 
  
  base_size = 30
  
  p = ggplot(plot.df, aes(x = variable, y = parameters)) + 
    geom_tile(aes(fill = value), colour = "white") + 
    scale_fill_gradientn(colours = c("#388697",'#F0F0F0',"#EB5160"),
                         na.value = "white", 
                         breaks=c(min(plot.df$value, na.rm = T), 0.5, max(plot.df$value, na.rm = T)),
                         labels=c('Higher value correlated with smaller\ndifferences in resistant carriers\nbetween the long and short duration ward',
                                  '\nNo effect\n', 
                                  'Higher value correlated with larger\ndifferences in resistant carriers\nbetween the long and short duration ward')) +
    theme_minimal(base_size = base_size) + 
    labs(x = "", y = "", fill = "") + 
    scale_y_discrete(expand = c(0, 0),
                     labels = c('','', '', 'Resistance transmission',
                                '','Antibiotic killing', 
                                '','','Within-host bacterial\ngrowth and decolonisation',
                                '','','Baseline carriage status', 
                                '','', '', 'Antibiotic prescriptions')) +
    scale_x_discrete(labels = c('', '', '', rep(c('Exclusive colonisation', 'Co-colonisation ', 'Within-host growth'), 6))) +
    theme(legend.position ='bottom', 
          legend.justification = c(0.6, 1),
          legend.key.height = unit(2,'cm'),
          axis.ticks = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(face='bold.italic', size = 12),
          axis.text.x = element_text(size = 12, angle = 330, hjust = 0, colour = "grey50"))+
    annotate(geom = 'text', x = 2, y = c(1:length(unique(plot.df$parameters))), 
             label = c('Bacterial interference factor', 'Amount of R transmitted', 
                       'R threshold to become R carrier', 'Transmission rate', 
                       'Killing of R', 'Killing of S', 
                       'Decolonisation rate',
                       'Growth rate of R', 'Growth rate of S', 
                       'Carrying capacity', '% Capacity occupied by GNB\non admission', 
                       '% R carriers on admission', 
                       '% Broad-spectrum AB prescribed\nduring admission*', 'Time to repeated AB prescription\nduring admission',
                       '% Broad-spectrum AB prescribed\nat admission*', '% Prescribed AB at admission*'), 
             colour ='grey40', size = 3.2) +
    guides(fill = guide_colorbar(nbin = 200, raster = F)) +
    geom_hline(yintercept=c(4.5, 6.5, 9.5, 12.5), color='grey20', size=0.5)+
    geom_vline(xintercept = c(1.5, 7.5, 13.5, 19.5)+2,  color = "black", size=0.75) + 
    geom_vline(xintercept = c(4.5, 10.5, 16.5)+2,  color = "grey", size=0.75) + 
    theme(plot.margin = unit(c(4, 2, 0, 0), "cm"), 
          text = element_text(size = 18), 
          legend.key.height = unit(0.75, "cm"),
          legend.key.width = unit(4, "cm"))
  
  outcomes = c('Resistance carriers amongst the\ntreated patients\n(Individual effect)',
               'Resistance carriers in the ward\n(Population effect)\n', 
               'Non-resistance carriers acquiring\nresistance carriage during admission\n(Population effect)')
  
  for (i in 1:length(outcomes))  {
    p = p + annotation_custom(
      grob = grid::textGrob(label = outcomes[i], hjust = 0.5, gp = gpar(cex = 1.2, fontface="bold")),
      ymin = 19.25,                   # Vertical position of the textGrob
      ymax = 20.25,
      xmin = c(6, 12, 18)[i],         # Note: The grobs are positioned outside the plot area
      xmax = c(7, 13, 19)[i])
  }
  
  eff.lab = 'Antibiotic active\nagainst susceptible\nand resistant organisms'
  ineff.lab = 'Antibiotic active\nonly against\nsusceptible organisms'
  
  abxr.outcomes = rep(c( ineff.lab, 
                         eff.lab), 3)
  x = 4 + 3 * c(0:5)
  for (i in 1:length(abxr.outcomes))  {
    p = p + annotation_custom(
      grob = grid::textGrob(label = abxr.outcomes[i], hjust = 0.5, gp = gpar(cex = 0.85)),
      ymin = 16.6,      # Vertical position of the textGrob
      ymax = 18.6,
      xmin = x[i],         # Note: The grobs are positioned outside the plot area
      xmax = (x+2)[i] )
  }
  
  
  # Code to override clipping
  gt = ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  grid.draw(gt)
  
}






