###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
################## show distribution of model outputs  ############################
###################################################################################

# combine data and output of each iteration 
make.df <- function(lhs.filename, outcome.type) {
  
  d = get(load(lhs.filename))  
  
  if(length(grep('treated', lhs.filename)) == 0) { # get res.names
    source('models/get_output_absdiff_simple3state.R')
  } else {
    source('models/get_output_treated_simple3state.R')
  }
  
  outcome.type.label = which(res.names == outcome.type)
  
  y = d$res[, outcome.type.label, 1]
  df = data.frame(model.name = as.character(d$call)[2], y = y)
  df$outcome.type = outcome.type

  df$model.name = as.factor(df$model.name)
  df$model.name = factor(df$model.name, 
                         levels = c('modelRun.simple', 'modelRun.cocarriage', 'modelRun.populationgrowth'), 
                         labels = c('Simple 3-state', 'Co-carriage 5-state', 'Population growth'))

  df$model.abxr.type = ifelse(max(d$data$abx.r) < 0.1, 'zero', 'notzero')
  df$model.abxr.type = as.factor(df$model.abxr.type)
  df$model.abxr.type = factor(df$model.abxr.type, levels = c('zero', 'notzero'))
  
  df$model.scenario = paste0(as.character(df$model.name), as.character(df$model.abxr.type), sep = ' ')
  
  return(df)
  
}

make.plot.d <- function(d.list){do.call('rbind', d.list)}

