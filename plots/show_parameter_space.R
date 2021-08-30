###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
###############################Check parameter space###############################
###################relationship between the parameters and outputs#################
###################################################################################
rm(list = ls(all.names = TRUE)) #clear environment 

source('plot_functions/plot_parameter_space.R')

### Ward characteristics
# how max.los affects los - max.los is generated from uniform distribution, 
# p.infect increasing from 0 to 1 - more prescriptions on admission - increases los 
# p.r.day1 increasing from 0 to 1 - more abx.r on admission, does not increase los 
(los.plot.max.los = plot.los(max.los.min = 3, max.los.max = 20, cum.r.1.min = 30, cum.r.1.max = 300, col = 'max.los'))
(los.plot.cum = plot.los(max.los.min = 3, max.los.max = 20, cum.r.1.min = 30, cum.r.1.max = 300, col = 'cum.r.1'))
(max.los.plot = plot.los(max.los.min = 3, max.los.max = 20, cum.r.1.min = 30, cum.r.1.max = 30, col = 'max.los'))
(cum.r.1.los.plot = plot.los(max.los.min = 7, max.los.max = 7, cum.r.1.min = 30, cum.r.1.max = 300, col = 'cum.r.1'))
ggarrange(los.plot.max.los, los.plot.cum, max.los.plot, cum.r.1.los.plot)
ggsave(filename = '../../../../Desktop/los.plot.jpeg', width = 40, height = 40, units = 'cm')

# how cum.r.1 affects use of broad spectrum antibiotics
(cum.r.1.hai.plot = cum.r.1.plot(cum.r.1.min=30, cum.r.1.max=500))
ggsave(filename = '../../../../Desktop/cum.r.1.hai.plot.jpeg', width = 15, height=10, units = 'cm')

### Baseline carriage status 
(baseline.simple = plot.baseline(model='simple'))
(baseline.binary = plot.baseline(model='binary'))
(baseline.frequency = plot.baseline(model='frequency'))
baseline.plot = ggarrange(baseline.simple, baseline.binary,baseline.frequency, nrow = 3)
ggsave(filename = '../../../../Desktop/baseline.plot.jpeg', width = 60, height = 50, units = 'cm')
### things to note here: 
# model 2: only prop_R corelates with S and R, the rest are hard to interpret anyways so only show prop_R 
# model 3: amount of max s, Sr, sr are adjusted to match frequency model! 

### Within host dynamics 
# model 2 repop.s 
(repop.s.plot = plot.repop (min=0.02, max=0.12, los=300, repop.name='repop.s', intercept = 45))
ggsave(filename = '../../../../Desktop/repop.s.plot.jpeg', width = 50, height=20, units = 'cm')

# model 2 repop.r
(repop.r.plot = plot.repop (min=0.02, max=0.2, repop.name='repop.r', intercept = 14, los = 30))
ggsave(filename = '../../../../Desktop/repop.r.plot.jpeg', width = 50, height=20, units = 'cm')

#transmission
(pi_ssr.plot = plot.pi_ssr(min=0.001, max=0.3, los=300))
ggsave(filename = '../../../../Desktop/pi_ssr.plot.jpeg', width = 50, height=20, units = 'cm')

# model 3 r_growth and s_growth 
min_s=0.03
max_s=0.3
min_r=0.3
max_r=1.4
growth(min_r=min_r, max_r=max_r, min_s=min_s, max_s=max_s, los=300, paratoexplore='prop_R')
growth(min_r=min_r, max_r=max_r, min_s=min_s, max_s=max_s, los=300, paratoexplore='r_thres')
growth(min_r=min_r, max_r=max_r, min_s=min_s, max_s=max_s, los=300, paratoexplore='total_prop')
growth(min_r=min_r, max_r=max_r, min_s=min_s, max_s=max_s, los=300, paratoexplore='K')
s.growth.plot=growth(min_r=min_r, max_r=max_r, min_s=min_s, max_s=max_s, los=300, paratoexplore='s_growth')
ggsave(filename = '../../../../Desktop/s.growth.plot.jpeg', width = 40, height=50, units = 'cm')
r.growth.plot=growth(min_r=min_r, max_r=max_r, min_s=min_s, max_s=max_s, los=300, paratoexplore='r_growth')
ggsave(filename = '../../../../Desktop/r.growth.plot.jpeg', width = 40, height=50, units = 'cm')

###Transmission and decolonisation 
#decolonisation
mu_plot = plot.mu(min=0.002, max=0.02, los=300)
ggsave(filename = '../../../../Desktop/mu_plot.jpeg', width = 50, height=20, units = 'cm')

#r_trans
# https://www.sciencedirect.com/science/article/pii/S0195670109000462 - transfer 0.15% for E. coli 
### hands contaminated to 7 × 108 cfu/mL, after 10min became 10000CFU
### after hand hygiene - 512500 cfu (range: 5750–3302500 cfu)
### median recovery from recipient higher from gloved hands (median 7500 cfu; range: 0–370000) 
####                                           bare hands (median 1000 cfu; range: 0–245000)
#https://www.ncbi.nlm.nih.gov/pubmed/26278471 - Enterobacteriaceae 1.2% of all microbiome of hands 
#https://www.ncbi.nlm.nih.gov/books/NBK144001/ - hands of HCWs 3.9 × 10^4 to 4.6 × 10^6 CFU/cm2

4*10^4*0.01*0.01 #exp(2)
4*10^6*0.015*0.02 #exp(7)

###Abx killing 
# models 1 and 2 
(abx.kill12.plot = abx.kill12 (min=0.1, max=0.5, los=14, abx.name='abx.r'))
ggsave(filename = '../../../../Desktop/abx.kill12.plot.jpeg', width = 50, height=20, units = 'cm')

#model 3
(abx.kill3.plot = abx.kill3(min_r=min_r, max_r=max_r, min_s=min_s, max_s=max_s, los=14, abx.name='abx.s', 
                            abx_min= 0.4, abx_max=0.8, paratoexplore='abx'))
ggsave(filename = '../../../../Desktop/abxs.kill3.plot.jpeg', width = 40, height=50, units = 'cm')
(abx.kill3.plot = abx.kill3(min_r = min_r, max_r = max_r, min_s = min_s, max_s = max_s, los = 30, abx.name='abx.r', 
                           abx_min = 1.3, abx_max= 1.7, paratoexplore='abx'))
ggsave(filename = '../../../../Desktop/abxr.kill3.plot.jpeg', width = 40, height=50, units = 'cm')

######################################################################################
###########Archive
#Shaw ISME 2019- microbiome recovery table 3 and 4 
d=data.frame(antibiotic=c('ciprofloxacin', 'clindamycin','amoxicillin','minocycline'), 
             Dmin=c(5.47,6.23,0.13,1.54), 
             Dmax=c(9.75,9.84,6.56,7.82),
             Amin=c(0.28,0.29,-0.66,-1.09),
             Amax=c(1.34,1.42,0.56,0.23), 
             phi1min=c(-0.69,-0.46,-1.96,-1.44),
             phi1max=c(0.16,0.34,0.31,1.29),
             phi2min=c(0.05,0.23,-1.58,1.01),
             phi2max=c(0.92,1.11,1.83,1.97))
t=1:100
D=1.34  #values in paper
A=-.03
phi1=-1.53 
phi2=0.09
#equation 3, figure 1c
y=((D*exp(phi1)*exp(phi2))/(exp(phi2)-exp(phi1)))*(exp(-exp(phi1*t))-exp(-exp(phi2*t)))+A*(1-exp(-exp(phi1*t)))
plot(t,y) #similar distribution as logistic growth

####K####
d=read.csv('gutdata/Ini_CTXm_copies_qPCR.csv') #SATURN gut data
mean(d$ini_16S_log) #total capacity in log, K
sd(d$ini_16S_log)
hist(d$ini_16S_log)

#####Gibson 2018: doubling time of e coli 15h 
# G (generation time) = (time, in minutes or hours)/n(number of generations), G = t/n
# GROWTH RATE = ln2/doubling time 
# log10(2)/(15/24) #15 hours doubling time - a daily growth rate of 0.5
# log10(2)/(25/24) #25 hours doubling time - lower limit a daily growth rate of 0.29 
# log10(2)/(5/24)  #5 hours doubling time - upper limit a daily growth rate of 1.44
# 
# t=1:500
# r_growth=0.29
# K=12
# r=200
# s=1000
# r_growth*r*(1 - (r + s)/exp(K)) #growth in 1 day 