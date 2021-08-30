########Show scatter of different parameters#############


scaleFUN <- function(x) sprintf("%.2f", x) #keep x axis decimal 2 places

scatter <-  function(lhs.data){
  
  # get the parameter values 
  sample = lhs.data$data
  
  # get out results 
  y.totalR = lhs.data$res[, 3, 1]
  y.newR = lhs.data$res[, 6, 1]
  
  # put them into dataframe 
  plotdata = data.frame(sample, y.totalR = y.totalR, y.newR = y.newR)
  
  colnames(plotdata)[which(colnames(plotdata) == 'n.bed')] = 'n'
  colnames(plotdata)[which(colnames(plotdata) == 'max.los')] = 'l'
  colnames(plotdata)[which(colnames(plotdata) == 'prop_R')] = 'p_R'
  colnames(plotdata)[which(colnames(plotdata) == 'prop_S')] = 'p_S'
  colnames(plotdata)[which(colnames(plotdata) == 'bif')] = 'b'
  colnames(plotdata)[which(colnames(plotdata) == 'pi_ssr')] = 'phi_s'
  colnames(plotdata)[which(colnames(plotdata) == 'repop.s')] = 'g_s'
  colnames(plotdata)[which(colnames(plotdata) == 'abx.s')] = 'alpha_s'
  colnames(plotdata)[which(colnames(plotdata) == 'abx.r')] = 'alpha_r'
  colnames(plotdata)[which(colnames(plotdata) == 'p.infect')] = 'omega_day1'
  colnames(plotdata)[which(colnames(plotdata) == 'p.r.day1')] = 'omega_day1.r'
  colnames(plotdata)[which(colnames(plotdata) == 'cum.r.1')] = 'omega_after'
  colnames(plotdata)[which(colnames(plotdata) == 'p.r.after')] = 'omega_after.r'
  colnames(plotdata)[which(colnames(plotdata) == 'long_dur')] = 't_long'
  colnames(plotdata)[which(colnames(plotdata) == 'short_dur')] = 't_short'
  
  colnames(plotdata)[which(colnames(plotdata) == 'prop_r')] = 'p_r'
  colnames(plotdata)[which(colnames(plotdata) == 'prop_Sr')] = 'p_Sr'
  colnames(plotdata)[which(colnames(plotdata) == 'fitness.r')] = 'f'
  
  colnames(plotdata)[which(colnames(plotdata) == 'r_thres')] = 'Gamma'
  colnames(plotdata)[which(colnames(plotdata) == 'r_trans')] = 'tau'
  colnames(plotdata)[which(colnames(plotdata) == 'total_prop')] = 'rho_e'
  colnames(plotdata)[which(colnames(plotdata) == 's_growth')] = 'c_s'
  
  plotdata.long.x = reshape2::melt(plotdata, # expand each sample variable 
                                   id.var = c('y.totalR', 'y.newR'),
                                   variable.name = 'x.variable',
                                   value.name = "x.value") 
  plotdata.long.x.y = reshape2::melt(plotdata.long.x, # expand each output variable 
                                     id.var = c('x.variable', 'x.value'), 
                                     variable.name = 'y.variable',
                                     value.name = "y.value")

  # plot 
  ggplot(plotdata.long.x.y, aes(x = x.value, y = y.value, group = y.variable, color = y.variable)) +
    geom_point(alpha = 0.2, size = 0.1) +
    geom_smooth(method = 'loess', formula = y ~ x, se = F) +
    facet_wrap(~ x.variable, scales = "free") +
    scale_x_continuous(labels = scaleFUN)+
    scale_color_manual(values = c('blue', 'orange'), 
                       name = 'Model output',
                       labels = c('Resistance carriers per bed per day', 
                                  'New Resistance acquisitions per admission')) + 
    labs(y = 'Model output', x = '') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = 'bottom')
}

