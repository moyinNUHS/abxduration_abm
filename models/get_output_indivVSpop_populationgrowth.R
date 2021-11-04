# get individual vs population effects

## create patient matrix where 1 column is one patient

source('models/los_abx_matrix_varydur.R')
source('models/summary_los.R')
source('models/model_populationgrowth.R')

amtperday.count <- function(df, capacity, abx.matrix, patient.matrix, max.dur){
  
  record = which(abx.matrix > 0)
  ## find who received abx (do not have to consider burn in)
  abx.pst = patient.matrix[record]
  
  amtperday.no = exp(df)/exp(capacity)
  amtperday.treated = amtperday.no[record]
  amtperday.treated.split = split(amtperday.treated, abx.pst)
  amtperday.treated.split.filled = lapply(amtperday.treated.split, `length<-`, max.dur + 1)
  
  amtperday.treated.matrix = do.call('cbind', amtperday.treated.split.filled)
  
  rowMeans(amtperday.treated.matrix, na.rm = T)
}

run_indivVSpop_populationgrowth <- function(p.infect, max.los, 
                                            K, total_prop, prop_R, pi_ssr, 
                                            r_trans, fitness.r,
                                            p.infect.after, p.r.day1, p.r.after,
                                            r_thres, s_growth, 
                                            abx.s, abx.r, iterations, max.dur){
  
  burn.in = 150
  n.day = 300
  n.bed = 50
  timestep = 1
  
  # empty matrix to store output - each row is a day, each col is a iteration 
  # store within host output - exact number of bacteria per day with each additional day of abx treatment 
  amtperday = list(eff = list(s = matrix(NA, nrow = max.dur+1, ncol = iterations),
                              r = matrix(NA, nrow = max.dur+1, ncol = iterations)), 
                   ineff = list(s = matrix(NA, nrow = max.dur+1, ncol = iterations),
                                r = matrix(NA, nrow = max.dur+1, ncol = iterations)))
  # store population output - number of R from within host selection and transmission with varying duration of abx applied to the ward 
  Rtype_iter = list(eff = array(NA, dim = c((n.day-burn.in), 4, iterations)), 
                    ineff = array(NA, dim = c((n.day-burn.in), 4, iterations)))
  Rtype_dur = list(eff = array(NA, dim = c((n.day-burn.in), 4,  max.dur)), 
                   ineff = array(NA, dim = c((n.day-burn.in), 4, max.dur)))
  
  
  for (dur in 0:max.dur){ # for each duration
    
    for(iter in 1:iterations){ # for each iteration 
      
      
      # get matrix of length of stay, abx prescribed, patients admitted
      matrixes = los_abx_table_varydur(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                                       p.infect = p.infect, p.r.day1 = p.r.day1, p.r.after = p.r.after, 
                                       p.infect.after = p.infect.after,  
                                       meanDur = dur, timestep=timestep)
      patient.matrix = matrixes[[1]]
      abx.matrix = matrixes[[2]]
      los.array = summary_los(patient.matrix = patient.matrix)
      
      # starting state for all the patients admitted 
      colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, total_prop=total_prop, prop_R=prop_R, r_thres=r_thres, K=K)
      
      # update values for each day - list 1 contains S, list 2 contains R
      colo.matrix_filled_iter = list()
      colo.matrix_filled_iter[['eff']] = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                                 pi_ssr=pi_ssr, total_prop=total_prop, K=K, 
                                                 fitness.r=fitness.r, r_thres=r_thres,r_trans=r_trans, s_growth=s_growth,
                                                 abx.s=abx.s, abx.r=abx.r, timestep=timestep)
      
      # update values for each day - list 1 contains S, list 2 contains R - if ineffective antibiotics given 
      colo.matrix_filled_iter[['ineff']] = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                                   pi_ssr=pi_ssr, total_prop=total_prop, K=K, 
                                                   fitness.r=fitness.r, r_thres=r_thres,r_trans=r_trans, s_growth=s_growth,
                                                   abx.s=abx.s, abx.r=0.00001, timestep=timestep)
      
      # summary of the output - timestep by bed in absolute numbers
      
      for (type in c('eff', 'ineff')){
        
        df.R = colo.matrix_filled_iter[[type]][[2]] # get matrix containing R 
        df.S = colo.matrix_filled_iter[[type]][[1]] # get matrix containing S 
        capacity = colo.matrix_filled_iter[[type]][[3]]
        df.S[which(df.S == -Inf)] = 0
        df.R[which(df.R == -Inf)] = 0
        
        if (dur == max.dur) {
          ##### WITHIN HOST
          
          # get vector of exact number of bacteria over `n.day` of abx treatment averaged over all treated patients
          amtperday[[type]][['s']][,iter] = amtperday.count(df.S, capacity = capacity, abx.matrix = abx.matrix, patient.matrix = patient.matrix, max.dur = max.dur) 
          amtperday[[type]][['r']][,iter] = amtperday.count(df.R, capacity = capacity, abx.matrix = abx.matrix, patient.matrix = patient.matrix, max.dur = max.dur)
        }
        
        ##### POPULATION 
        r_thres_matrix = colo.matrix_filled_iter[[type]][[4]]
        colo_vector_filled_iter_TF = df.R >= r_thres_matrix # turn values into T (R carrier) or F (did not reach threshold for R carrier)
        df.trans = colo.matrix_filled_iter[[type]][[5]] # transmission event table 
        df.Rtype = matrix(NA, nrow = nrow(colo_vector_filled_iter_TF), ncol = ncol(colo_vector_filled_iter_TF)) # fill in all NA as stand-ins
        
        ## type 1 - admitted as R 
        df.Rtype[which(colo.matrix[[2]] >= r_thres_matrix)] = 'Admitted as R'
        ## type 2 - admitted as S
        df.Rtype[which(colo.matrix[[2]] < r_thres_matrix)] = 'S'
        
        # classify other types of R
        if (any(colo_vector_filled_iter_TF == T)) { # if there are any Rs
          
          for (d in 2:nrow(colo_vector_filled_iter_TF)){# for each row 
            
            tofill = is.na(df.Rtype[d,])
            
            ## classify each R into transmitted, or within host selection case
            R = colo_vector_filled_iter_TF[d, ] == T
            oldR = colo_vector_filled_iter_TF[d-1, ] == T
            newR = R & !oldR
            
            # if R is new, transmitted or within host 
            if (any(newR) == T) { # if there are new Rs 
              trans = df.trans[d, ] == T
              df.Rtype[d, trans & newR & tofill] = 'Transmitted'
              
              select = df.trans[d, ] != T
              df.Rtype[d, select & newR & tofill] = 'Selection'
            }
            
            # if R the previous day, follow the previous day's reason
            df.Rtype[d, which(R & oldR & tofill)] = df.Rtype[d-1, which(R & oldR &tofill)]
            
            # if not filled, then S 
            df.Rtype[d, is.na(df.Rtype[d,])] = 'S'
          }
        }
        
        # proportion of R per day 
        Rtype.dur.iter.type = t(apply(df.Rtype, 1, function(x){ #sum of every row for each type
          c('Admitted as R' = sum(x == 'Admitted as R', na.rm = T), 
            'Transmitted' = sum(x == 'Transmitted', na.rm = T), 
            'Selection' = sum(x == 'Selection', na.rm = T), 
            'S' = sum(x == 'S', na.rm = T))
        }))
        
        # Discard first `burn_in` days as burn-in
        Rtype_iter[[type]][, , iter] = Rtype.dur.iter.type[(burn.in+1):n.day,]
        
      }
      
    }
    
    for (type in c('eff', 'ineff')){ # get mean across iterations
      Rtype_dur[[type]][, , dur] =  apply(Rtype_iter[[type]], c(1,2), mean)
    }
    
  }
  
  
  ##########################
  ## clean data for plot 
  ##########################
  
  #### WITHIN HOST 
  amtperday$eff$s = rowMeans(amtperday$eff$s)
  amtperday$eff$r = rowMeans(amtperday$eff$r)
  amtperday$ineff$s = rowMeans(amtperday$ineff$s)
  amtperday$ineff$r = rowMeans(amtperday$ineff$r)
  
  amtperday.wide = do.call('cbind.data.frame', amtperday)
  amtperday.wide$day = 0:max.dur
  amtperday.long = reshape2::melt(amtperday.wide, id.var = 'day')
  amtperday.long$abx = gsub("\\..*","",amtperday.long$variable)
  amtperday.long$res = gsub("*..","",amtperday.long$variable)
  
  #### POPULATION
  Rtype.eff.mean = t(apply(Rtype_dur$eff, c(2,3), mean, na.rm = T)) #mean across each ward in each observational period
  Rtype.eff.mean = as.data.frame(Rtype.eff.mean)
  Rtype.eff.mean$day = 0: (max.dur -1)
  colnames(Rtype.eff.mean) = c('Admitted as R', 'Transmitted', 'Selection', 'Non-carrier', 'Day')
  Rtype.eff = reshape2::melt(Rtype.eff.mean, id.var = 'Day')
  Rtype.eff$abx = 'Effective antibiotic available'
  
  Rtype.ineff.mean = t(apply(Rtype_dur$ineff, c(2,3), mean, na.rm = T))
  Rtype.ineff.mean = as.data.frame(Rtype.ineff.mean)
  Rtype.ineff.mean$day = 0: (max.dur -1)
  colnames(Rtype.ineff.mean) = c('Admitted as R', 'Transmitted', 'Selection', 'Non-carrier', 'Day')
  Rtype.ineff = reshape2::melt(Rtype.ineff.mean, id.var = 'Day')
  Rtype.ineff$abx = 'Effective antibiotic not available'
  
  Rtype.long = rbind.data.frame(Rtype.eff, Rtype.ineff)
  
  return(list(withinhost = amtperday.long, pop = Rtype.long))
  
}

