###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
#### Check monotonicity of the parameters - Hoeffding's D and Spearman's ##########
###################################################################################
scaleFUN <- function(x) sprintf("%.2f", x) #keep x axis decimal 2 places

get.hd.spm <- function(x, y){
  
  hd = scaleFUN(hoeffd(x[['sample']], x[[y]])$D[1,2])
  hd.p = scaleFUN(hoeffd(x[['sample']], x[[y]])$P[1,2])
  
  spm = scaleFUN(cor.test(x[['sample']], x[[y]], method = "spearman")$estimate)
  spm.p = scaleFUN(cor.test(x[['sample']], x[[y]], method = "spearman")$p.value)
  
  return(c(hd = hd, hd.p = hd.p, spm = spm, spm.p = spm.p))
}


# produce table 
monotonicity.tab <- function(lhs.data) {
  
  # get the parameter values 
  sample = lhs.data$data
  
  # get out results 
  y.totalR = lhs.data$res[, 3, 1]
  y.newR = lhs.data$res[, 6, 1]
  
  # put them into list
  lhs.data.list = list()
  for (i in 1:ncol(sample)){
    lhs.data.list[[colnames(sample[i])]] = list(sample = sample[ ,i], 
                                                y.totalR = y.totalR, 
                                                y.newR = y.newR)
  }
  
  # calculate Hoeffding's D and Spearman's
  out = lapply(lhs.data.list, function(x){
    
    calculated.values = c(
      
      ## total R 
      get.hd.spm(x, 'y.totalR'),
      
      ## new R 
      get.hd.spm(x, 'y.newR')
    )
    
    raw = matrix(calculated.values, nrow = 2, byrow = T)
    raw[grep('0.00', raw)] = '0'
    
    raw[,2][which(raw[,2] == '0')] = '<0.01'
    raw[,4][which(raw[,4] == '0')] = '<0.01'
    
    return(raw)
    
  })
  
  # present as table 
  out.df = do.call('rbind.data.frame', out)
  out.df = cbind.data.frame(rep(colnames(sample), each = 2), # parameter names 
                            rep(c('R carriers per day', 'New R acquisitions per admission'), ncol(sample)), 
                            out.df)
  
  return(out.df)
}

get.latex <- function(df, model){
  
  rle.lengths <- rle(df[[1]])$lengths
  first <- !duplicated(df[[1]])
  df[[1]][!first] <- ""
  
  # define appearance of \multirow
  df[[1]][first] = paste0("\\midrule\\multirow{", rle.lengths, "}{*}{\\textbf{", df[[1]][first], "}}")
  
  # caption
  strCaption <- paste0("Hoeffding's D measure and Spearman's rank correlation measures the ", model)
  
  # set up xtable output
  print(xtable(df, digits = c(0, 0, 0, 3, 1, 0, 6), # first zero "represents" row numbers which we skip later
               align = "lllrr|rr",  # align and put a vertical line (first "l" again represents column of row numbers)
               caption = strCaption, label = paste('monotonicity', model)),
        size = "footnotesize", #Change size; useful for bigger tables "normalsize" "footnotesize"
        include.rownames = FALSE, #Don't print rownames
        include.colnames = FALSE, #We create them ourselves
        caption.placement = "top", #"top", NULL
        hline.after = NULL, #We don't need hline; we use booktabs
        floating = TRUE, # whether \begin{Table} should be created (TRUE) or not (FALSE)
        sanitize.text.function = force, # Important to treat content of first column as latex function
        add.to.row = list(pos = list(-1, 2,
                                     nrow(df)),
                          command = c(paste("\\toprule \n",  # NEW row
                                            "\\multicolumn{2}{c}{} & \\multicolumn{2}{c}{\\textbf{Hoeffding`s D measure}} & \\multicolumn{2}{c}{\\textbf{Spearman`s rank correlation measure}} \\\\\n",
                                            "\\cmidrule(l){3-4} \\cmidrule(l){5-6}\n",
                                            " & & Measure & p-value & Measure & p-value \\\\\n", # NEW row 
                                            "\\midrule \n"
                          ),
                          paste("\\cmidrule(l){3-4} \\cmidrule(l){5-6}\n" # we may also use 'pos' and 'command' to add a midrule
                          ),
                          paste("\\bottomrule \n"  # paste is used as it is more flexible regarding adding lines
                          )
                          )
        )
  )
  
}

get.csv <- function(df.simple, df.cocarriage, df.populationgrowth){
  
  df.list = list(df.simple, df.cocarriage, df.populationgrowth)
  
  clean.list = lapply(df.list, function(df){
    colnames(df) = c('parameters', 'outcome.type', 'HD', 'HD.p', 'SPM', 'SPM.p')
    df$HD = paste0(df$HD, ' (', df$HD.p, ')')
    df$SPM = paste0(df$SPM, ' (', df$SPM.p, ')')
    
    df$parameters[which(df$parameter == 'n.bed')] = 'n'
    df$parameters[which(df$parameter == 'max.los')] = 'l'
    df$parameters[which(df$parameter == 'prop_R')] = 'p_R'
    df$parameters[which(df$parameter == 'prop_r')] = 'p_r'
    df$parameters[which(df$parameter == 'prop_S')] = 'p_S'
    df$parameters[which(df$parameter == 'prop_Sr')] = 'p_Sr'
    df$parameters[which(df$parameter == 'bif')] = 'b'
    
    df$parameters[which(df$parameter == 'pi_ssr')] = 'phi_s'
    df$parameters[which(df$parameter == 'repop.s')] = 'g_s'
    df$parameters[which(df$parameter == 'abx.s')] = 'alpha_s'
    df$parameters[which(df$parameter == 'abx.r')] = 'alpha_r'
    df$parameters[which(df$parameter == 'p.infect')] = 'omega_day1'
    df$parameters[which(df$parameter == 'p.r.day1')] = 'omega_day1.r'
    df$parameters[which(df$parameter == 'cum.r.1')] = 'omega_after'
    df$parameters[which(df$parameter == 'p.r.after')] = 'omega_after.r'
    df$parameters[which(df$parameter == 'long_dur')] = 't_long'
    df$parameters[which(df$parameter == 'short_dur')] = 't_short'
    
    df$parameters[which(df$parameter == 'fitness.r')] = 'f'
    
    df$parameters[which(df$parameter == 'r_thres')] = 'Gamma'
    df$parameters[which(df$parameter == 'r_trans')] = 'tau'
    df$parameters[which(df$parameter == 'total_prop')] = 'rho_e'
    df$parameters[which(df$parameter == 's_growth')] = 'c_s'
    
    data.frame(Parameter = df$parameters, 
               Outcome = df$outcome.type,
               HD = df$HD, 
               SPM = df$SPM)
  })

  out = clean.list %>% reduce(full_join, by = c("Parameter", "Outcome"))
  out[is.na(out)] = '-'
  
  out$Parameter[duplicated(out$Parameter)] = ''
  colnames(out)[3:ncol(out)] = rep(c('HD', 'SPM'), 3)

  return(out)

}
