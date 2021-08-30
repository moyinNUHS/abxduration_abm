# functions to determine number of iterations per run

get_AA <- function(model, abx.r.upper, iterationstotry) {
  
  # SOURCE MODELS AND SAMPLE PARAMETER SPACE 
  if (model == 'simple3state'){
    
    source('models/get_output_absdiff_simple3state.R')
    
    parameters <- c(
      n.bed = c(runif(1, min=5, max=50)),         # "n.bed", number of beds in the ward
      max.los = c(runif(1, min=3, max=20)),       # "max.los", mean of length of stay
      prop_R = c(runif(1, min=0.01, max=0.8)),    # "prop_R",probability of initial carriage of resistant organisms
      prop_S = c(runif(1, min=0, max=1)),         # "prop_S", proportion of S in the population of S and ss
      bif = c(runif(1, min=0, max=1)),            # "bif", bacterial interference factor
      pi_ssr = c(runif(1, min=0.001, max=0.3)),   # "pi_ssr" probability of being transmitted r to ss (ss—> ssr)
      repop.s = c(runif(1, min=0.02, max=0.12)),  # "repop.s1" probability of ss repopulated to S (Palleja, Nature Biology, 2018 on gut recovery ~9 months)
      mu = c(runif(1, min=0.002, max=0.02)),      # "mu, probability of decolonisation (Haggai Bar-Yoseph, JAC, 2016, decreasing colonization rates from 76.7% (95% CI=69.3%–82.8%) at 1 month to 35.2% (95% CI=28.2%–42.9%) at 12 months of follow-up)
      abx.s = c(runif(1, min=0.1, max=0.5)),      # "abx.s", probability of S becoming ss after being on narrow spectrum antibiotics
      abx.r =  c(runif(1, min=0, max=0.5)),       # "abx.r", probability of R becoming ss after being on broad spectrum antibiotics
      p.infect = c(runif(1, min=0.1, max=1)),     # "p.infect", probability of being prescribed antibiotics
      cum.r.1 = c(runif(1, min=30, max= 300)),    # admission day when cummulative probability of HAI requiring abx.r is 1
      p.r.day1 = c(runif(1, min=0.1, max=1)),     # probability of being prescribed broad spectrum antibiotic on admission 
      p.r.after = c(runif(1, min=0.1, max=1)),    #  probability of being prescribed broad spectrum antibiotic after admission
      short_dur = c(runif(1, min=3, max=7)),      # "short_dur", mean short duration of antibiotics (normal distribution)
      long_dur = c(runif(1, min=14, max=21))      # "long_dur", mean long duration of antibiotics (normal distribution)
    )
    
  } else if (model == 'cocarriage5state') {
    
    source('models/get_output_absdiff_cocarriage5state.R')
    
    parameters <- c(
      n.bed = c(runif(1,min=5, max=50)),         #n.bed; number of beds in the ward
      max.los = c(runif(1,min=3, max=20)),       #max.los; mean of length of stay (exponential distribution)
      prop_R = c(runif(1,min=0, max=0.8)),       #probability of initial carriage of resistant organisms
      prop_r = c(runif(1,min=0, max=1)),         #proportion of S in (S+s): prob_start_S <- prop_S_nonR*(1-prob_R)
      prop_Sr = c(runif(1,min=0, max=1)),        #proportion of Sr in (r+R): prob_start_Sr <- prop_Sr_inR*prob_R
      prop_S = c(runif(1,min=0, max=1)),         #proportion of sr in (r+r): prob_start_sr <- prop_sr_inR*prob_R
      bif = c(runif(1,min=0, max=1)),            #bacterial interference factor (pi_ssr = pi_r1 * bif )
      pi_ssr = c(runif(1,min=0.001, max=0.3)),   #probability of being transmitted r to ss (ss—> ssr)
      repop.s = c(runif(1,min=0.02, max=0.12)),  #probability of regrowth of S  (s—>S)
      fitness.cost = c(runif(1,min=0, max=1)),  #probability of regrowth of s (sr—> sR)
      mu = c(runif(1,min=0.002, max=0.02)),      #probability of being decolonised to S (Sr—> S) 
      abx.s = c(runif(1,min=0.1, max=0.5)),      #probability of clearing S to become s
      abx.r = c(runif(1,min=0.1, max=0.5)),      #probability of clearing R to become r
      p.infect = c(runif(1,min=0.1, max=1)),     #probability of being prescribed narrow spectrum antibiotic
      cum.r.1 = c(runif(1,min=30, max=300)),     #admission day when cumulative prabability of HAI requiring abx.r is 1
      p.r.day1 = c(runif(1,min=0.1, max=1)),     #probability of being prescribed broad spectrum antibiotic on day 1 of admission 
      p.r.after = c(runif(1, min=0.1, max=1)),    #  probability of being prescribed broad spectrum antibiotic after admission
      short_dur = c(runif(1,min=3, max=7)),      #mean short duration of antibiotics (normal distribution) 
      long_dur = c(runif(1,min=14, max=21))      #mean long duration of antibiotics (normal distribution) 
    )
    
  } else if (model == 'populationgrowth') {
    
    source('models/get_output_absdiff_populationgrowth.R')
    
    parameters <- c(
      n.bed = c(runif(1,min=5, max=50)),         #n.bed; number of beds in the ward
      max.los = c(runif(1,min=3, max=20)),       #max.los; mean of length of stay (normal distribution)
      p.infect = c(runif(1,min=0.1, max=1)),     #probability of being prescribed narrow spectrum antibiotic
      cum.r.1 = c(runif(1,min=30, max=300)),     #admission day when cummulative prabability of HAI requiring abx.r is 1
      p.r.day1 = c(runif(1,min=0.1, max=1)),     #probability of being prescribed broad spectrum antibiotic on day 1 of admission 
      p.r.after = c(runif(1, min=0.1, max=1)),    #  probability of being prescribed broad spectrum antibiotic after admission
      K = c(runif(1,min=18, max=24)),            # gut holding capacity, on log scale, largest R number possible is exp(300) - typical colonic bacteria 10^14 number/mL content https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4991899/
      total_prop = c(runif(1,min=0.1, max=0.9)), # mean of total starting amount of e coli on log scale
      prop_R = c(runif(1,min=0,max=0.8)),        # mean of starting amount of resistant gut bacteria on log scale
      pi_ssr = c(runif(1,min=0.001,max=0.3)),    # pi_ssr = daily probability of transmitting resistant E coli
      r_trans = c(runif(1,min=0.01,max=0.2)),         # r_mean= R threshold level for tranmissibility
      fitness.cost = c(runif(1,min=0.3, max=1.4)),   # fitness.cost = growth constant for logistic growth
      r_thres = c(runif(1,min=0.01, max=0.2)),       # r_thres = threshold amount of bacteria before R can be transmitted
      s_growth = c(runif(1,min=0.03,max=0.3)),   # fitness.cost = growth constant for logistic growth
      abx.s = c(runif(1,min=0.1,max=0.8)),       # abxr_killr = amount of r killed by broad spectrum abx r
      abx.r = c(runif(1,min=0.1,max=0.8)),       # abxr_kills = amount of s killed by broad spectrum abx r
      short_dur = c(runif(1,min=3, max=7)),      #mean short duration of narrow spectrum antibiotics (normal distribution) 
      long_dur = c(runif(1,min=14, max=21))      #mean long duration of narrow spectrum antibiotics (normal distribution)
    )
    
  }
  if (abx.r.upper < 0.1) {parameters[['abx.r']] = 0.00000000001}
  
  # number of runs to perform
  # iterationstotry = iterations we are going to test 
  numberofrepeatsineachiteration = 20
  no.runs = max(iterationstotry) * numberofrepeatsineachiteration
  
  # Consistency analysis operates by contrasting distributions of simulation responses, all generated using the same fixed set of parameter
  # values and containing identical numbers of simulation samples.
  params.matrix = matrix(rep(parameters, no.runs), nrow = no.runs, byrow = T)
  colnames(params.matrix) = names(parameters)
  
  # run models 
  if (model == 'simple3state'){
    aa_data = t(apply(params.matrix, 1, function(values) {
      run_absdiff_simple3state(n.bed = values[['n.bed']], max.los = values[['max.los']], 
                               prop_R = values[['prop_R']], prop_S = values[['prop_S']], 
                               bif = values[['bif']], pi_ssr = values[['pi_ssr']], repop.s = values[['repop.s']], mu = values[['mu']], 
                               abx.s = values[['abx.s']], abx.r = values[['abx.r']], p.infect = values[['p.infect']], 
                               cum.r.1 = values[['cum.r.1']], p.r.day1 = values[['p.r.day1']], p.r.after = values[['p.r.after']], short_dur = values[['short_dur']], long_dur = values[['long_dur']])
    }))
  } else if (model == 'cocarriage5state') {
    aa_data = t(apply(params.matrix, 1, function(values) {
      run_absdiff_cocarriage5state(n.bed = values[['n.bed']], max.los = values[['max.los']], 
                                   prop_R = values[['prop_R']], prop_r = values[['prop_r']], prop_Sr = values[['prop_Sr']], prop_S = values[['prop_S']],
                                   bif = values[['bif']], pi_ssr = values[['pi_ssr']], repop.s = values[['repop.s']], fitness.cost = values[['fitness.cost']], mu = values[['mu']], 
                                   abx.s = values[['abx.s']], abx.r = values[['abx.r']], p.infect = values[['p.infect']], 
                                   cum.r.1 = values[['cum.r.1']], p.r.day1 = values[['p.r.day1']], p.r.after = values[['p.r.after']], short_dur = values[['short_dur']], long_dur = values[['long_dur']])
    }))
    
  } else if (model == 'populationgrowth') {
    aa_data = t(apply(params.matrix, 1, function(values) {
      run_absdiff_populationgrowth(n.bed = values[['n.bed']], max.los = values[['max.los']], 
                                   p.infect = values[['p.infect']], cum.r.1 = values[['cum.r.1']],
                                   p.r.day1 = values[['p.r.day1']], p.r.after = values[['p.r.after']], K = values[['K']], total_prop = values[['total_prop']], 
                                   prop_R = values[['prop_R']], pi_ssr = values[['pi_ssr']], r_trans = values[['r_trans']], fitness.cost = values[['fitness.cost']], 
                                   r_thres = values[['r_thres']], s_growth = values[['s_growth']], abx.s = values[['abx.s']], abx.r = values[['abx.r']],
                                   short_dur = values[['short_dur']], long_dur = values[['long_dur']])
    }))
  }
  out = cbind(params.matrix, aa_data)
  colnames(out) = c(colnames(params.matrix), res.names)
  
  #store simulation results in appropriate folders - one output per folder 
  abx.label = ifelse(abx.r.upper < 0.1, 'abxrZERO', 'abxrNOTZERO')
  dirtostoreAAruns = paste0(working.directory.main, 'AA_tests/AA_runs/', model, '/', abx.label) # make sure this directory is already created! 
  
  for (i in iterationstotry){
    
    data = as.data.frame(out) # data to store
    setwd(dirtostoreAAruns)    # go to the folder to store the data 
    numbertoputinfolderi = i * numberofrepeatsineachiteration
    dir.create(as.character(i)) # create folders for `iterationstotry`
    
    for (k in 1 : numberofrepeatsineachiteration){
      
      setwd(paste0(dirtostoreAAruns, '/', as.character(i))) # go to the folders for the `iterationstotry`
      dir.create(as.character(k)) # create a folder for the repeat in the iteration 
      
      for (g in 1:i){
        setwd(paste0(dirtostoreAAruns, '/', as.character(i), '/', as.character(k)))
        dir.create(as.character(g))
        setwd(paste0(dirtostoreAAruns, '/', as.character(i), '/', as.character(k), '/', as.character(g)))
        write.csv(data[1,], 'aa_data_simple.csv', row.names = FALSE)
        data = data[-1,]
        
      }
    }
    setwd(working.directory.main)
  }
  
  #Running Aleatory Analysis 
  FILEPATH <- dirtostoreAAruns #already in dirtostoreAAruns stated above 
  # Sample sizes (number of simulation replicates in each distribution) to be analysed
  SAMPLESIZES <- iterationstotry
  # The simulation output measures to be analysed
  MEASURES <- res.names
  # Number of distributions being compared. Default: 20, as performed by Read et al
  NUMSUBSETSPERSAMPLESIZE <- numberofrepeatsineachiteration
  # Output file name containing the simulation responses.
  RESULTFILENAME <- "aa_data_simple.csv"
  # Notes the column in the CSV results file where the results start.
  OUTPUTFILECOLSTART <- ncol(params.matrix) + 1
  # Last column of the output measure results
  OUTPUTFILECOLEND <- ncol(out) 
  # The A-Test value either side of 0.5 which should be considered a 'large difference'
  # between two sets of results. Use of 0.23 was taken from the Vargha-Delaney publication
  LARGEDIFFINDICATOR <- 0.23
  # A-Test values above 0.5 (no difference) which should be considered as small,
  # medium, and large differences between two result sets. Used in the graph
  # summarising all sample sizes.
  SMALL <- 0.56
  MEDIUM <- 0.66
  LARGE <- 0.73
  
  # A summary file is created containing the median A-Test values for each sample size.
  SUMMARYFILENAME <- "AA_ATestMaxAndMedians.csv"
  
  #calculate median values for each sample size
  aa_summariseReplicateRuns(FILEPATH = FILEPATH, SAMPLESIZES = SAMPLESIZES, MEASURES = MEASURES, 
                            RESULTFILENAME = RESULTFILENAME, # Output file name containing the simulation responses
                            NUMSUBSETSPERSAMPLESIZE = NUMSUBSETSPERSAMPLESIZE,
                            OUTPUTFILECOLSTART = OUTPUTFILECOLSTART, OUTPUTFILECOLEND = OUTPUTFILECOLEND, 
                            SUMMARYFILENAME = SUMMARYFILENAME) #A summary file containing the median A-Test values for each sample size
  
  # The results of the A-Test comparisons of the twenty subsets for each sample size
  # are stored within an output file. 
  ATESTRESULTSFILENAME <- "AA_ATest_Scores.csv"
  
  # calculate A-test scores - get a csv file of the A test scores and graphs output for each sample size 
  a_test_results <- aa_getATestResults(FILEPATH = FILEPATH, SAMPLESIZES = SAMPLESIZES, NUMSUBSETSPERSAMPLESIZE = NUMSUBSETSPERSAMPLESIZE, 
                                       MEASURES = MEASURES, ATESTRESULTSFILENAME = ATESTRESULTSFILENAME, LARGEDIFFINDICATOR = LARGEDIFFINDICATOR, 
                                       AA_SIM_RESULTS_FILE = "AA_ATestMaxAndMedians.csv")
  
  # A dataframe of summary of the max and median A test scores for each sample size 
  sample_summary <- aa_sampleSizeSummary(FILEPATH, SAMPLESIZES, MEASURES, SUMMARYFILENAME, ATESTRESULTS_OBJECT = a_test_results)
  
  # Graphs, by the sample_summary object - gives a summary graph of the sample_summary object above 
  # Name of the graph which summarises the analysis results for all sample sizes.
  GRAPHOUTPUTFILE <- "AA_ATestMaxes.pdf"
  aa_graphSampleSizeSummary(FILEPATH, MEASURES, 300, SMALL, MEDIUM, LARGE, GRAPHOUTPUTFILE, SAMPLESUMMARY_OBJECT = sample_summary)
  
}