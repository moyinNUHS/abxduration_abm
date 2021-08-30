###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
################## Plot partial correlation coefficient ###########################
###################################################################################

source('plots/plot_functions/plot_prcc.R')

simple_0=get(load('runs/LHSdiff_simple_370_zero23Apr2020_2000BST.Rdata'))
binary_0=get(load('runs/LHSdiff_binary_370_zero_24Apr2020_0012BST.Rdata'))
freq_0=get(load('runs/LHSdiff_frequency_370_zero_23Apr2020_2243BST.Rdata'))

simple_not0=get(load('runs/LHSdiff_simple_370_notzero23Apr2020_2041BST.Rdata'))
binary_not0=get(load('runs/LHSdiff_binary_370_notzero_23Apr2020_2328BST.Rdata'))
freq_not0=get(load('runs/LHSdiff_frequency_370_notzero_23Apr2020_2143BST.Rdata'))

simple0 = plotprcc(simple_0)
binary0 = plotprcc(binary_0)
freq0 = plotprcc(freq_0)

ggarrange(simple0, binary0, freq0, nrow=3, labels = c('Model 1', 'Model 2', 'Model 3'))

simplenot0 = plotprcc(LHS.simple)
binarynot0 = plotprcc(binary_not0)
freqnot0 = plotprcc(freq_not0)

ggarrange(simplenot0, binarynot0, freqnot0, nrow=3, labels = c('Model 1', 'Model 2', 'Model 3'))
