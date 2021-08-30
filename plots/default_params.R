############## Default parameters ##################

### FOR STABILITY GRAPH 
#common parameters 
common.para = list(n.bed = 20,       # n.bed= number of beds in the ward
                   max.los = 20,     # mean.max.los= mean of max length of stay (exponential distribution)
                   n.day = 400,      # always 400 for stability graph
                   short_dur = 5,
                   long_dur = 15,
                   prop_R = 0.4,     # Probability of being colonized with resistant strain on admission
                   pi_ssr = 0.03,    # pi_ssr= probability of R transmitting 
                   cum.r.1 = 100, 
                   p.infect = 0.5,   # p=probability of receiving antibiotic
                   p.r.day1 = 0.5, 
                   p.r.after = 0.5,
                   timestep = 1)


specific.para = list(
  ############# Simple 3-state model #################
  simple3state =  list(prop_S = 0.3,     # Proportion of large S within non-resistant states (S+s)
                       bif = 0.8,         # bacterial interference factor 
                       mu = 0.01,         # mu= probability of clearance of Sr to become S
                       repop.s = 0.06, 
                       abx.s = 0.3,
                       abx.r = 0.3), 
  
  ############## Cocarriage 5-state model ###################
  cocarriage5state = list(prop_S = 0.3,     # Proportion of large S within non-resistant states (S+s)
                          prop_Sr = 0.4,                  
                          prop_r = 0.5,    
                          bif = 0.8,         # bacterial interference factor
                          mu = 0.01,         # mu= probability of clearance of Sr to become S
                          repop.s = 0.06, 
                          fitness.r = 0.5,     # probability of repopulation of Sr to become sR 
                          abx.s = 0.3,
                          abx.r = 0.3),
  
  ############# Population growth model #################
  populationgrowth = list(K = exp(22),           # gut holding capacity
                          total_prop = 0.9,  # mean of total starting amount of gut bacteria on log scale
                          fitness.r = 0.5,    # fitness.r = growth cost for R
                          s_growth = 0.2,
                          r_trans = 0.1,       # r_thres = R threshold level for tranmissibility
                          r_thres = 0.1, 
                          abx.s = 0.7,
                          abx.r = 0.7)
)


