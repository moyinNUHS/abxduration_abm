substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# get ranking position of each parameter, add NA if a parameter is not in a particular model
getposition <- function(lhs.filename, labs = labs.df, outcome.type){
  
  d = get(load(lhs.filename))  
  
  if(length(grep('treated', lhs.filename)) == 0) { # get res.names
    source('models/get_output_absdiff_simple3state.R')
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
  
  names(out) = rev(param.labels)
  
  return(out)
}

# clean data for plot 
clean_ranks_data <- function(treated, newR, totalR){
  
  
  # combine the data 
  d = c(treated, totalR, newR)
  ranks.df.wide = as.data.frame(do.call('cbind', d))
  models = c('simple', 'cocarriage', 'populationgrowth')
  colnames(ranks.df.wide) = c(paste0(models, '_notzero_', 'treated'), 
                              paste0(models, '_zero_', 'treated'),
                              paste0(models, '_notzero_', 'totalR'), 
                              paste0(models, '_zero_', 'totalR'),
                              paste0(models, '_notzero_', 'newR'), 
                              paste0(models, '_zero_', 'newR'))
  ranks.df.wide$parameters = rownames(ranks.df.wide)
  ranks.df.wide$empty = NA # create empty first column to place parameter names 
  ranks.df = reshape2::melt(ranks.df.wide, id.var = c('parameters'))
  ranks.df$variable = as.factor(ranks.df$variable)
  ranks.df$variable = factor(ranks.df$variable, 
                             levels = c('empty', 
                                        paste0(models, '_zero_', 'treated'), 
                                        paste0(models, '_zero_', 'totalR'),
                                        paste0(models, '_zero_', 'newR'), 
                                        paste0(models, '_notzero_', 'treated'),
                                        paste0(models, '_notzero_', 'totalR'), 
                                        paste0(models, '_notzero_', 'newR')))
  ranks.df$value[intersect(which(ranks.df$parameters == 'alpha_r'), grep('_zero', ranks.df$variable))] = NA
  ranks.df$parameters = as.factor(ranks.df$parameters)
  ranks.df$parameters = factor(ranks.df$parameters, 
                               levels = rev(param.labels))
  
  #remove parameters that need not be shown 
  to.remove = c('t_short', 't_long', 'p_S', 'p_Sr', 'p_r')
  plot.df = ranks.df[-which(ranks.df$parameters %in% to.remove),]
  
  return(plot.df)
}


# plot function 
plot_paraheatmap <- function (plot.df) { 
  
  base_size = 9
  
  p = ggplot(plot.df, aes(x = variable, y = parameters)) + 
    geom_tile(aes(fill = value), colour = "white") + 
    scale_fill_gradientn(colours = c("#388697",'#fbfae5',"#EB5160"),
                         na.value = "white", 
                         breaks=c(min(plot.df$value, na.rm = T), 0.5, max(plot.df$value, na.rm = T)),
                         labels=c('Higher value decreases resistant\ncarriers in the long duration ward',
                                  'No effect', 
                                  'Higher value increases resistant\ncarriers in the long duration ward')) +
    theme_minimal(base_size = base_size) + 
    labs(x = "", y = "", fill = "") + 
    scale_y_discrete(expand = c(0, 0),
                     labels = c('','Antibiotic killing', 
                                '','','','Bacterial growth\nand decolonisation',
                                '','','Baseline carriage status', 
                                '','', '', 'Resistance transmission',
                                '','', '', 'Antibiotic prescriptions',
                                '', 'Ward characteristics')) +
    scale_x_discrete(labels = c('', rep(c('Simple 3-state', 'Co-carriage 5-state', 'Population growth'), 6))) +
    theme(legend.position ='bottom', 
          legend.justification = c(0.6, 1),
          legend.key.width = unit(2,'cm'),
          axis.ticks = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(face='bold.italic', size = 12),
          axis.text.x = element_text(size = 12, angle = 330, hjust = 0, colour = "grey50"))+
    annotate(geom = 'text', x = 1, y = c(1:length(unique(plot.df$parameters))), 
             label = unique(plot.df$parameters), 
             colour ='grey40', size = 3) +
    guides(fill = guide_colorbar(nbin = 200, raster = F))+
    geom_hline(yintercept=c(2.5, 6.5, 9.5, 13.5, 17.5), color='grey20', size=0.5)+
    geom_vline(xintercept = c(1.5, 7.5, 13.5, 19.5),  color = "black", size=0.75) + 
    geom_vline(xintercept = c(4.5, 10.5, 16.5),  color = "grey", size=0.75) + 
    theme(plot.margin = unit(c(4, 2, 0, 0), "cm"), 
          text = element_text(size = 18), 
          legend.key.width = unit(5, "cm"))
  
  outcomes = c('Proportion of resistance carriers\namongst the treated patients\n',
               'Difference in proportion of\nresistance carriers per day\n', 
               'Difference in proportion of\nnon-resistance carriers acquiring\nresistance carriage during admission')
  
  for (i in 1:length(outcomes))  {
    p = p + annotation_custom(
      grob = grid::textGrob(label = outcomes[i], hjust = 0.5, gp = gpar(cex = 1.2, fontface="bold")),
      ymin = 22,      # Vertical position of the textGrob
      ymax = 23,
      xmin = c(4, 10, 16)[i],         # Note: The grobs are positioned outside the plot area
      xmax = c(5, 11, 17)[i])
  }
  
  abxr.outcomes = rep(c('Treated with\nineffective antibiotic', 
                        'Treated with\neffective antibiotic'), 3)
  x = 2 + 3 * c(0:5)
  for (i in 1:length(abxr.outcomes))  {
    p = p + annotation_custom(
      grob = grid::textGrob(label = abxr.outcomes[i], hjust = 0.5, gp = gpar(cex = 1)),
      ymin = 19.5,      # Vertical position of the textGrob
      ymax = 21.5,
      xmin = x[i],         # Note: The grobs are positioned outside the plot area
      xmax = (x+2)[i] )
  }
  
  
  # Code to override clipping
  gt = ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  grid.draw(gt)
  
}






