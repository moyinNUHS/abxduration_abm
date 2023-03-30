###########################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###########
####### run models for magnitude of change in resistance prevalence #######################
###########################################################################################
rm(list = ls()) # Clean working environment

source('magnitude_allmodels.R')

#run data 

cre.abxr = c(abx.r.min = 0, abx.r.max = 0.0001)
high.trans.prev = c(pi_ssr.min = 0.3, pi_ssr.max= 0.3001, prop_R.min = 0.8, prop_R.max = 0.8001, p.r.day1.min = 0.8, p.r.day1.max = 0.8001)
low.trans.prev = c(pi_ssr.min = 0.001, pi_ssr.max= 0.0011, prop_R.min = 0.1, prop_R.max = 0.1001, p.r.day1.min = 0.1, p.r.day1.max = 0.1001)
high.prescript = c(p.infect.min = 0.8, p.infect.max = 0.8001, cum.r.1.min = 30, cum.r.1.max = 100)
low.prescript = c(p.infect.min = 0.1, p.infect.max = 0.1001, cum.r.1.min = 200, cum.r.1.max = 300)

### ### ### 
### CRE ###
### ### ### 

### Scenario 1 high transmission/ prevalence, low prescription 
run_magnitude_simple(pi_ssr.min = high.trans.prev['pi_ssr.min'], pi_ssr.max = high.trans.prev['pi_ssr.max'], 
                     prop_R.min = high.trans.prev['prop_R.min'], prop_R.max = high.trans.prev['prop_R.max'], 
                     p.r.day1.min = high.trans.prev['p.r.day1.min'], p.r.day1.max = high.trans.prev['p.r.day1.max'],
                     abx.r.min = cre.abxr['abx.r.min'], abx.r.max = cre.abxr['abx.r.max'], 
                     abx.s.min = 0.5, abx.s.max = 0.5001, 
                     p.infect.min = low.prescript['p.infect.min'], p.infect.max = low.prescript['p.infect.max'], 
                     cum.r.1.min = low.prescript['cum.r.1.min'], cum.r.1.max = low.prescript['cum.r.1.max'])

run_magnitude_binary(pi_ssr.min = high.trans.prev['pi_ssr.min'], pi_ssr.max = high.trans.prev['pi_ssr.max'], 
                     prop_R.min = high.trans.prev['prop_R.min'], prop_R.max = high.trans.prev['prop_R.max'], 
                     p.r.day1.min = high.trans.prev['p.r.day1.min'], p.r.day1.max = high.trans.prev['p.r.day1.max'],
                     abx.r.min = cre.abxr['abx.r.min'], abx.r.max = cre.abxr['abx.r.max'], 
                     abx.s.min = 0.5, abx.s.max = 0.5001, 
                     p.infect.min = low.prescript['p.infect.min'], p.infect.max = low.prescript['p.infect.max'], 
                     cum.r.1.min = low.prescript['cum.r.1.min'], cum.r.1.max = low.prescript['cum.r.1.max'])

run_magnitude_frequency(pi_ssr.min = high.trans.prev['pi_ssr.min'], pi_ssr.max = high.trans.prev['pi_ssr.max'], 
                        prop_R.min = high.trans.prev['prop_R.min'], prop_R.max = high.trans.prev['prop_R.max'], 
                        p.r.day1.min = high.trans.prev['p.r.day1.min'], p.r.day1.max = high.trans.prev['p.r.day1.max'],
                        abx.r.min = cre.abxr['abx.r.min'], abx.r.max = cre.abxr['abx.r.max'], 
                        abx.s.min = 0.8, abx.s.max = 0.8001, 
                        p.infect.min = low.prescript['p.infect.min'], p.infect.max = low.prescript['p.infect.max'], 
                        cum.r.1.min = low.prescript['cum.r.1.min'], cum.r.1.max = low.prescript['cum.r.1.max'])

### Scenario 2 high transmission/ prevalence, high prescription 
run_magnitude_simple(pi_ssr.min = high.trans.prev['pi_ssr.min'], pi_ssr.max = high.trans.prev['pi_ssr.max'], 
                     prop_R.min = high.trans.prev['prop_R.min'], prop_R.max = high.trans.prev['prop_R.max'], 
                     p.r.day1.min = high.trans.prev['p.r.day1.min'], p.r.day1.max = high.trans.prev['p.r.day1.max'],
                     abx.r.min = cre.abxr['abx.r.min'], abx.r.max = cre.abxr['abx.r.max'], 
                     abx.s.min = 0.5, abx.s.max = 0.5001, 
                     p.infect.min = high.prescript['p.infect.min'], p.infect.max = high.prescript['p.infect.max'], 
                     cum.r.1.min = high.prescript['cum.r.1.min'], cum.r.1.max = high.prescript['cum.r.1.max'])

run_magnitude_binary(pi_ssr.min = high.trans.prev['pi_ssr.min'], pi_ssr.max = high.trans.prev['pi_ssr.max'], 
                     prop_R.min = high.trans.prev['prop_R.min'], prop_R.max = high.trans.prev['prop_R.max'], 
                     p.r.day1.min = high.trans.prev['p.r.day1.min'], p.r.day1.max = high.trans.prev['p.r.day1.max'],
                     abx.r.min = cre.abxr['abx.r.min'], abx.r.max = cre.abxr['abx.r.max'], 
                     abx.s.min = 0.5, abx.s.max = 0.5001, 
                     p.infect.min = high.prescript['p.infect.min'], p.infect.max = high.prescript['p.infect.max'], 
                     cum.r.1.min = high.prescript['cum.r.1.min'], cum.r.1.max = high.prescript['cum.r.1.max'])

run_magnitude_frequency(pi_ssr.min = high.trans.prev['pi_ssr.min'], pi_ssr.max = high.trans.prev['pi_ssr.max'], 
                        prop_R.min = high.trans.prev['prop_R.min'], prop_R.max = high.trans.prev['prop_R.max'], 
                        p.r.day1.min = high.trans.prev['p.r.day1.min'], p.r.day1.max = high.trans.prev['p.r.day1.max'],
                        abx.r.min = cre.abxr['abx.r.min'], abx.r.max = cre.abxr['abx.r.max'], 
                        abx.s.min = 0.8, abx.s.max = 0.8001, 
                        p.infect.min = high.prescript['p.infect.min'], p.infect.max = high.prescript['p.infect.max'], 
                        cum.r.1.min = high.prescript['cum.r.1.min'], cum.r.1.max = high.prescript['cum.r.1.max'])

### Scenario 3 low transmission/ prevalence, low prescription 
run_magnitude_simple(pi_ssr.min = low.trans.prev['pi_ssr.min'], pi_ssr.max = low.trans.prev['pi_ssr.max'], 
                     prop_R.min = low.trans.prev['prop_R.min'], prop_R.max = low.trans.prev['prop_R.max'], 
                     p.r.day1.min = low.trans.prev['p.r.day1.min'], p.r.day1.max = low.trans.prev['p.r.day1.max'],
                     abx.r.min = cre.abxr['abx.r.min'], abx.r.max = cre.abxr['abx.r.max'], 
                     abx.s.min = 0.5, abx.s.max = 0.5001, 
                     p.infect.min = low.prescript['p.infect.min'], p.infect.max = low.prescript['p.infect.max'], 
                     cum.r.1.min = low.prescript['cum.r.1.min'], cum.r.1.max = low.prescript['cum.r.1.max'])

run_magnitude_binary(pi_ssr.min = low.trans.prev['pi_ssr.min'], pi_ssr.max = low.trans.prev['pi_ssr.max'], 
                     prop_R.min = low.trans.prev['prop_R.min'], prop_R.max = low.trans.prev['prop_R.max'], 
                     p.r.day1.min = low.trans.prev['p.r.day1.min'], p.r.day1.max = low.trans.prev['p.r.day1.max'],
                     abx.r.min = cre.abxr['abx.r.min'], abx.r.max = cre.abxr['abx.r.max'], 
                     abx.s.min = 0.5, abx.s.max = 0.5001, 
                     p.infect.min = low.prescript['p.infect.min'], p.infect.max = low.prescript['p.infect.max'], 
                     cum.r.1.min = low.prescript['cum.r.1.min'], cum.r.1.max = low.prescript['cum.r.1.max'])

run_magnitude_frequency(pi_ssr.min = low.trans.prev['pi_ssr.min'], pi_ssr.max = low.trans.prev['pi_ssr.max'], 
                        prop_R.min = low.trans.prev['prop_R.min'], prop_R.max = low.trans.prev['prop_R.max'], 
                        p.r.day1.min = low.trans.prev['p.r.day1.min'], p.r.day1.max = low.trans.prev['p.r.day1.max'],
                        abx.r.min = cre.abxr['abx.r.min'], abx.r.max = cre.abxr['abx.r.max'], 
                        abx.s.min = 0.8, abx.s.max = 0.8001, 
                        p.infect.min = low.prescript['p.infect.min'], p.infect.max = low.prescript['p.infect.max'], 
                        cum.r.1.min = low.prescript['cum.r.1.min'], cum.r.1.max = low.prescript['cum.r.1.max'])

### Scenario 4 low transmission/ prevalence, high prescription 
run_magnitude_simple(pi_ssr.min = low.trans.prev['pi_ssr.min'], pi_ssr.max = low.trans.prev['pi_ssr.max'], 
                     prop_R.min = low.trans.prev['prop_R.min'], prop_R.max = low.trans.prev['prop_R.max'], 
                     p.r.day1.min = low.trans.prev['p.r.day1.min'], p.r.day1.max = low.trans.prev['p.r.day1.max'],
                     abx.r.min = cre.abxr['abx.r.min'], abx.r.max = cre.abxr['abx.r.max'], 
                     abx.s.min = 0.5, abx.s.max = 0.5001, 
                     p.infect.min = high.prescript['p.infect.min'], p.infect.max = high.prescript['p.infect.max'], 
                     cum.r.1.min = high.prescript['cum.r.1.min'], cum.r.1.max = high.prescript['cum.r.1.max'])

run_magnitude_binary(pi_ssr.min = low.trans.prev['pi_ssr.min'], pi_ssr.max = low.trans.prev['pi_ssr.max'], 
                     prop_R.min = low.trans.prev['prop_R.min'], prop_R.max = low.trans.prev['prop_R.max'], 
                     p.r.day1.min = low.trans.prev['p.r.day1.min'], p.r.day1.max = low.trans.prev['p.r.day1.max'],
                     abx.r.min = cre.abxr['abx.r.min'], abx.r.max = cre.abxr['abx.r.max'], 
                     abx.s.min = 0.5, abx.s.max = 0.5001, 
                     p.infect.min = high.prescript['p.infect.min'], p.infect.max = high.prescript['p.infect.max'], 
                     cum.r.1.min = high.prescript['cum.r.1.min'], cum.r.1.max = high.prescript['cum.r.1.max'])

run_magnitude_frequency(pi_ssr.min = low.trans.prev['pi_ssr.min'], pi_ssr.max = low.trans.prev['pi_ssr.max'], 
                        prop_R.min = low.trans.prev['prop_R.min'], prop_R.max = low.trans.prev['prop_R.max'], 
                        p.r.day1.min = low.trans.prev['p.r.day1.min'], p.r.day1.max = low.trans.prev['p.r.day1.max'],
                        abx.r.min = cre.abxr['abx.r.min'], abx.r.max = cre.abxr['abx.r.max'], 
                        abx.s.min = 0.8, abx.s.max = 0.8001, 
                        p.infect.min = high.prescript['p.infect.min'], p.infect.max = high.prescript['p.infect.max'], 
                        cum.r.1.min = high.prescript['cum.r.1.min'], cum.r.1.max = high.prescript['cum.r.1.max'])

### ### ### ### 
###  3GCRE  ###
### ### ### ### 

### Scenario 1 high transmission/ prevalence, low prescription 
run_magnitude_simple(pi_ssr.min = high.trans.prev['pi_ssr.min'], pi_ssr.max = high.trans.prev['pi_ssr.max'], 
                     prop_R.min = high.trans.prev['prop_R.min'], prop_R.max = high.trans.prev['prop_R.max'], 
                     p.r.day1.min = high.trans.prev['p.r.day1.min'], p.r.day1.max = high.trans.prev['p.r.day1.max'],
                     abx.r.min = 0.5, abx.r.max = 0.5001, 
                     abx.s.min = 0.5, abx.s.max = 0.5001, 
                     p.infect.min = low.prescript['p.infect.min'], p.infect.max = low.prescript['p.infect.max'], 
                     cum.r.1.min = low.prescript['cum.r.1.min'], cum.r.1.max = low.prescript['cum.r.1.max'])

run_magnitude_binary(pi_ssr.min = high.trans.prev['pi_ssr.min'], pi_ssr.max = high.trans.prev['pi_ssr.max'], 
                     prop_R.min = high.trans.prev['prop_R.min'], prop_R.max = high.trans.prev['prop_R.max'], 
                     p.r.day1.min = high.trans.prev['p.r.day1.min'], p.r.day1.max = high.trans.prev['p.r.day1.max'],
                     abx.r.min = 0.5, abx.r.max = 0.5001, 
                     abx.s.min = 0.5, abx.s.max = 0.5001, 
                     p.infect.min = low.prescript['p.infect.min'], p.infect.max = low.prescript['p.infect.max'], 
                     cum.r.1.min = low.prescript['cum.r.1.min'], cum.r.1.max = low.prescript['cum.r.1.max'])

run_magnitude_frequency(pi_ssr.min = high.trans.prev['pi_ssr.min'], pi_ssr.max = high.trans.prev['pi_ssr.max'], 
                        prop_R.min = high.trans.prev['prop_R.min'], prop_R.max = high.trans.prev['prop_R.max'], 
                        p.r.day1.min = high.trans.prev['p.r.day1.min'], p.r.day1.max = high.trans.prev['p.r.day1.max'],
                        abx.r.min = 0.8, abx.r.max = 0.8001, 
                        abx.s.min = 0.8, abx.s.max = 0.8001, 
                        p.infect.min = low.prescript['p.infect.min'], p.infect.max = low.prescript['p.infect.max'], 
                        cum.r.1.min = low.prescript['cum.r.1.min'], cum.r.1.max = low.prescript['cum.r.1.max'])

### Scenario 2 high transmission/ prevalence, high prescription 
run_magnitude_simple(pi_ssr.min = high.trans.prev['pi_ssr.min'], pi_ssr.max = high.trans.prev['pi_ssr.max'], 
                     prop_R.min = high.trans.prev['prop_R.min'], prop_R.max = high.trans.prev['prop_R.max'], 
                     p.r.day1.min = high.trans.prev['p.r.day1.min'], p.r.day1.max = high.trans.prev['p.r.day1.max'],
                     abx.r.min = 0.5, abx.r.max = 0.5001, 
                     abx.s.min = 0.5, abx.s.max = 0.5001, 
                     p.infect.min = high.prescript['p.infect.min'], p.infect.max = high.prescript['p.infect.max'], 
                     cum.r.1.min = high.prescript['cum.r.1.min'], cum.r.1.max = high.prescript['cum.r.1.max'])

run_magnitude_binary(pi_ssr.min = high.trans.prev['pi_ssr.min'], pi_ssr.max = high.trans.prev['pi_ssr.max'], 
                     prop_R.min = high.trans.prev['prop_R.min'], prop_R.max = high.trans.prev['prop_R.max'], 
                     p.r.day1.min = high.trans.prev['p.r.day1.min'], p.r.day1.max = high.trans.prev['p.r.day1.max'],
                     abx.r.min = 0.5, abx.r.max = 0.5001, 
                     abx.s.min = 0.5, abx.s.max = 0.5001, 
                     p.infect.min = high.prescript['p.infect.min'], p.infect.max = high.prescript['p.infect.max'], 
                     cum.r.1.min = high.prescript['cum.r.1.min'], cum.r.1.max = high.prescript['cum.r.1.max'])

run_magnitude_frequency(pi_ssr.min = high.trans.prev['pi_ssr.min'], pi_ssr.max = high.trans.prev['pi_ssr.max'], 
                        prop_R.min = high.trans.prev['prop_R.min'], prop_R.max = high.trans.prev['prop_R.max'], 
                        p.r.day1.min = high.trans.prev['p.r.day1.min'], p.r.day1.max = high.trans.prev['p.r.day1.max'],
                        abx.r.min = 0.8, abx.r.max = 0.8001, 
                        abx.s.min = 0.8, abx.s.max = 0.8001, 
                        p.infect.min = high.prescript['p.infect.min'], p.infect.max = high.prescript['p.infect.max'], 
                        cum.r.1.min = high.prescript['cum.r.1.min'], cum.r.1.max = high.prescript['cum.r.1.max'])

### Scenario 3 low transmission/ prevalence, low prescription 
run_magnitude_simple(pi_ssr.min = low.trans.prev['pi_ssr.min'], pi_ssr.max = low.trans.prev['pi_ssr.max'], 
                     prop_R.min = low.trans.prev['prop_R.min'], prop_R.max = low.trans.prev['prop_R.max'], 
                     p.r.day1.min = low.trans.prev['p.r.day1.min'], p.r.day1.max = low.trans.prev['p.r.day1.max'],
                     abx.r.min = 0.5, abx.r.max = 0.5001, 
                     abx.s.min = 0.5, abx.s.max = 0.5001, 
                     p.infect.min = low.prescript['p.infect.min'], p.infect.max = low.prescript['p.infect.max'], 
                     cum.r.1.min = low.prescript['cum.r.1.min'], cum.r.1.max = low.prescript['cum.r.1.max'])

run_magnitude_binary(pi_ssr.min = low.trans.prev['pi_ssr.min'], pi_ssr.max = low.trans.prev['pi_ssr.max'], 
                     prop_R.min = low.trans.prev['prop_R.min'], prop_R.max = low.trans.prev['prop_R.max'], 
                     p.r.day1.min = low.trans.prev['p.r.day1.min'], p.r.day1.max = low.trans.prev['p.r.day1.max'],
                     abx.r.min = 0.5, abx.r.max = 0.5001, 
                     abx.s.min = 0.5, abx.s.max = 0.5001, 
                     p.infect.min = low.prescript['p.infect.min'], p.infect.max = low.prescript['p.infect.max'], 
                     cum.r.1.min = low.prescript['cum.r.1.min'], cum.r.1.max = low.prescript['cum.r.1.max'])

run_magnitude_frequency(pi_ssr.min = low.trans.prev['pi_ssr.min'], pi_ssr.max = low.trans.prev['pi_ssr.max'], 
                        prop_R.min = low.trans.prev['prop_R.min'], prop_R.max = low.trans.prev['prop_R.max'], 
                        p.r.day1.min = low.trans.prev['p.r.day1.min'], p.r.day1.max = low.trans.prev['p.r.day1.max'],
                        abx.r.min = 0.8, abx.r.max = 0.8001, 
                        abx.s.min = 0.8, abx.s.max = 0.8001, 
                        p.infect.min = low.prescript['p.infect.min'], p.infect.max = low.prescript['p.infect.max'], 
                        cum.r.1.min = low.prescript['cum.r.1.min'], cum.r.1.max = low.prescript['cum.r.1.max'])

### Scenario 4 low transmission/ prevalence, high prescription 
run_magnitude_simple(pi_ssr.min = low.trans.prev['pi_ssr.min'], pi_ssr.max = low.trans.prev['pi_ssr.max'], 
                     prop_R.min = low.trans.prev['prop_R.min'], prop_R.max = low.trans.prev['prop_R.max'], 
                     p.r.day1.min = low.trans.prev['p.r.day1.min'], p.r.day1.max = low.trans.prev['p.r.day1.max'],
                     abx.r.min = 0.5, abx.r.max = 0.5001, 
                     abx.s.min = 0.5, abx.s.max = 0.5001, 
                     p.infect.min = high.prescript['p.infect.min'], p.infect.max = high.prescript['p.infect.max'], 
                     cum.r.1.min = high.prescript['cum.r.1.min'], cum.r.1.max = high.prescript['cum.r.1.max'])

run_magnitude_binary(pi_ssr.min = low.trans.prev['pi_ssr.min'], pi_ssr.max = low.trans.prev['pi_ssr.max'], 
                     prop_R.min = low.trans.prev['prop_R.min'], prop_R.max = low.trans.prev['prop_R.max'], 
                     p.r.day1.min = low.trans.prev['p.r.day1.min'], p.r.day1.max = low.trans.prev['p.r.day1.max'],
                     abx.r.min = 0.5, abx.r.max = 0.5001, 
                     abx.s.min = 0.5, abx.s.max = 0.5001, 
                     p.infect.min = high.prescript['p.infect.min'], p.infect.max = high.prescript['p.infect.max'], 
                     cum.r.1.min = high.prescript['cum.r.1.min'], cum.r.1.max = high.prescript['cum.r.1.max'])

run_magnitude_frequency(pi_ssr.min = low.trans.prev['pi_ssr.min'], pi_ssr.max = low.trans.prev['pi_ssr.max'], 
                        prop_R.min = low.trans.prev['prop_R.min'], prop_R.max = low.trans.prev['prop_R.max'], 
                        p.r.day1.min = low.trans.prev['p.r.day1.min'], p.r.day1.max = low.trans.prev['p.r.day1.max'],
                        abx.r.min = 0.8, abx.r.max = 0.8001, 
                        abx.s.min = 0.8, abx.s.max = 0.8001, 
                        p.infect.min = high.prescript['p.infect.min'], p.infect.max = high.prescript['p.infect.max'], 
                        cum.r.1.min = high.prescript['cum.r.1.min'], cum.r.1.max = high.prescript['cum.r.1.max'])
