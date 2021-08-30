###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
################################## Plot cobweb ####################################
###################################################################################
source('plot_functions/plot_cobweb.R')

output=get(load('runs/LHSdiff_simple_355_notzero15Apr2020_2214BST.Rdata'))

jpeg("../../../../Desktop/cobweb.jpeg", width = 1000, height = 750)
plot_cobweb(output)
dev.off()
