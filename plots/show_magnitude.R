###########################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###########
############ show magnitude of change in resistance prevalence ############################
###########################################################################################
rm(list=ls()) # Clean working environment

source('plot_functions/plot_magnitude.R')

# load plot data 
########### CRE
rawlist.cre.hh = list(simple = get(load('runs/magnitude_simple_200_abxr0_pissr0.3_pinfect0.803Jun2020_1836BST.Rdata')),
                      binary = get(load('runs/magnitude_binary_200_abxr0_pissr0.3_pinfect0.803Jun2020_1847BST.Rdata')),
                      freq = get(load('runs/magnitude_frequency_200_abxr0_pissr0.3_pinfect0.803Jun2020_1902BST.Rdata')))
crehh = get_plotdata(rawlist = rawlist.cre.hh, 'crehh')

rawlist.cre.hl = list(simple = get(load('runs/magnitude_simple_200_abxr0_pissr0.3_pinfect0.103Jun2020_1759BST.Rdata')),
                      binary = get(load('runs/magnitude_binary_200_abxr0_pissr0.3_pinfect0.103Jun2020_1809BST.Rdata')),
                      freq = get(load('runs/magnitude_frequency_200_abxr0_pissr0.3_pinfect0.103Jun2020_1825BST.Rdata')))
crehl = get_plotdata(rawlist = rawlist.cre.hl, 'crehl')

rawlist.cre.lh = list(simple = get(load('runs/magnitude_simple_200_abxr0_pissr0.001_pinfect0.803Jun2020_1948BST.Rdata')),
                      binary = get(load('runs/magnitude_binary_200_abxr0_pissr0.001_pinfect0.803Jun2020_1958BST.Rdata')),
                      freq = get(load('runs/magnitude_frequency_200_abxr0_pissr0.001_pinfect0.803Jun2020_2016BST.Rdata')))
crelh = get_plotdata(rawlist = rawlist.cre.lh, 'crelh')

rawlist.cre.ll = list(simple = get(load('runs/magnitude_simple_200_abxr0_pissr0.001_pinfect0.103Jun2020_1911BST.Rdata')),
                      binary = get(load('runs/magnitude_binary_200_abxr0_pissr0.001_pinfect0.103Jun2020_1921BST.Rdata')),
                      freq = get(load('runs/magnitude_frequency_200_abxr0_pissr0.001_pinfect0.103Jun2020_1939BST.Rdata')))
crell = get_plotdata(rawlist = rawlist.cre.ll, 'crell')


########### 3GCRE

rawlist.3gcre.hh = list(simple = get(load('runs/magnitude_simple_200_abxr0.5_pissr0.3_pinfect0.803Jun2020_2059BST.Rdata')),
                      binary = get(load('runs/magnitude_binary_200_abxr0.5_pissr0.3_pinfect0.803Jun2020_2110BST.Rdata')),
                      freq = get(load('runs/magnitude_frequency_200_abxr0.8_pissr0.3_pinfect0.803Jun2020_2126BST.Rdata')))
tgcrehh = get_plotdata(rawlist = rawlist.3gcre.hh, '3gcrehh')

rawlist.3gcre.hl = list(simple = get(load('runs/magnitude_simple_200_abxr0.5_pissr0.3_pinfect0.103Jun2020_2024BST.Rdata')),
                      binary = get(load('runs/magnitude_binary_200_abxr0.5_pissr0.3_pinfect0.103Jun2020_2034BST.Rdata')),
                      freq = get(load('runs/magnitude_frequency_200_abxr0.8_pissr0.3_pinfect0.103Jun2020_2050BST.Rdata')))
tgcrehl = get_plotdata(rawlist = rawlist.3gcre.hl, '3gcrehl')

rawlist.3gcre.lh = list(simple = get(load('runs/magnitude_simple_200_abxr0.5_pissr0.001_pinfect0.803Jun2020_2213BST.Rdata')),
                      binary = get(load('runs/magnitude_binary_200_abxr0.5_pissr0.001_pinfect0.803Jun2020_2224BST.Rdata')),
                      freq = get(load('runs/magnitude_frequency_200_abxr0.8_pissr0.001_pinfect0.803Jun2020_2243BST.Rdata')))
tgcrelh = get_plotdata(rawlist = rawlist.3gcre.lh, '3gcrelh')

rawlist.3gcre.ll = list(simple = get(load('runs/magnitude_simple_200_abxr0.5_pissr0.001_pinfect0.103Jun2020_2135BST.Rdata')),
                      binary = get(load('runs/magnitude_binary_200_abxr0.5_pissr0.001_pinfect0.103Jun2020_2145BST.Rdata')),
                      freq = get(load('runs/magnitude_frequency_200_abxr0.8_pissr0.001_pinfect0.103Jun2020_2204BST.Rdata')))
tgcrell = get_plotdata(rawlist = rawlist.3gcre.ll, '3gcrell')

#list of plot data 
plotlist.cre = list(crehh, crehl, crelh, crell)
plotlist.3gcre = list(tgcrehh, tgcrehl, tgcrelh, tgcrell)

# plot graphs 
plot_magnitude(plotlist.left = plotlist.3gcre, plotlist.right = plotlist.cre, n.bed = 15, n.days = 30)
ggsave('/Users/moyin/Desktop/show_magnitude.pdf', units = 'cm', width = 30, height = 15)
