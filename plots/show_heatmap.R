###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
###### Plot heatmap exploration matrices of pairs of parameters ###################
###################################################################################

rm(list=ls()) # Clean working environment

###################### Adjust parameter ranges: START ##################
model = 'simple'
source("default_params.R")
source(paste0('model_', model, '.R'))
pixels = 15  #how many pixels per heatmap

params_ranges_simple = list(esbl = list(seq(0, 0.6, length.out = pixels), #repop.s
                                        seq(0, 0.05, length.out = pixels),#mu
                                        seq(0, 0.3, length.out = pixels), #pi_ssr
                                        seq(0, 1, length.out = pixels),   #bif
                                        seq(0, 0.3, length.out = pixels), #abx.s
                                        seq(0, 0.3, length.out = pixels)), #abx.r
                            cpe = list(seq(0, 0.6, length.out = pixels), #repop.s
                                       seq(0, 0.01, length.out = pixels),#mu
                                       seq(0, 0.1, length.out = pixels), #pi_ssr
                                       seq(0, 1, length.out = pixels), #bif
                                       seq(0, 0.3, length.out = pixels))) #abx.s 

params_ranges_binary = list(esbl = list(seq(0, 0.6, length.out = pixels), #repop.s
                                        seq(0, 0.1, length.out = pixels), #repop.r  
                                        seq(0, 0.03, length.out = pixels),#mu
                                        seq(0, 0.3, length.out = pixels), #pi_ssr
                                        seq(0, 1, length.out = pixels),   #bif
                                        seq(0, 0.3, length.out = pixels), #abx.s
                                        seq(0, 0.3, length.out = pixels)), #abx.r
                            cpe = list(seq(0, 0.6, length.out = pixels), #repop.s
                                       seq(0, 0.1, length.out = pixels), #repop.r  
                                       seq(0, 0.03, length.out = pixels),#mu
                                       seq(0, 0.3, length.out = pixels), #pi_ssr
                                       seq(0, 1, length.out = pixels),   #bif
                                       seq(0, 0.3, length.out = pixels)))#abx.s

params_ranges_frequency = list(esbl = list(seq(0, 0.5, length.out = pixels), #s_growth
                                           seq(0, 1.5, length.out = pixels), #r_growth
                                           seq(6, 10, length.out = pixels), #r_thres
                                           seq(0, 5, length.out = pixels), #r_trans
                                           seq(0, 0.3, length.out = pixels), #pi_ssr
                                           seq(0.5, 0.9, length.out = pixels), #abx.s
                                           seq(0.5, 0.9, length.out = pixels)), #abx.r
                               cpe = list(seq(0, 0.5, length.out = pixels), #s_growth
                                          seq(0, 1.5, length.out = pixels), #r_growth
                                          seq(6, 10, length.out = pixels), #r_thres
                                          seq(0, 5, length.out = pixels), #r_trans
                                          seq(0, 0.3, length.out = pixels), #pi_ssr
                                          seq(0.5, 0.9, length.out = pixels))) #abx.s

plots_data = data_heatmap(model = model, params_ranges = params_ranges_binary)
save(plots_data, file=paste0("./runs/heatmap2by2_", model,Sys.time(), ".Rdata"))

###################### Adjust parameter ranges: END ##################


# plot heatmap
source('plot_functions/plot_heatmap.R')
plots_data = get(load(paste0('./runs/heatmap2by2_simple.Rdata')))
plot_heatmap (model = 'simple', plots_data = plots_data)
ggsave(filename = paste0('../../../../../Desktop/heatmap_simple.pdf'), 
       width = 30, height = 20, units = 'cm')

plots_data = get(load(paste0('./runs/heatmap2by2_binary.Rdata')))
plot_heatmap (model = 'binary', plots_data = plots_data)
ggsave(filename = paste0('../../../../../Desktop/heatmap_binary.pdf'), 
       width = 30, height = 20, units = 'cm')


