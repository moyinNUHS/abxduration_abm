
####Model 1 and 2 



########repop.r: returning traveller 25% colonised with esbl after 1-6 weeks of travel
t=1:50
repopr=0.0005
#at time t, the probability of s-->S 
becoming_R_at_time_t= (1-repopr)^(t-1)*repopr
cumsum_becoming_R_at_time_t=cumsum(becoming_R_at_time_t)
plot(t, cumsum_becoming_R_at_time_t)
#repopr=0.05 for ~ 1 week for cummulative risk to become 0.25
#repopr=0.0005 for ~ 6 weeks for cummulative risk to become 0.25

########pi_ssr
t=1:40
pi_ssr=0.05
#at time t, the probability of s-->S 
becoming_R_at_time_t= (1-pi_ssr)^(t-1)*pi_ssr
cumsum_becoming_R_at_time_t=cumsum(becoming_R_at_time_t)
plot(t, cumsum_becoming_R_at_time_t)
#repopr=0.02 for ~ 50 days for cummulative risk to become 0.60 -hilty 2012 cid
#repopr=0.05 for  ~ 15 days for cummulative risk to become 0.50

########timestep
#prob_r competes with abx which goes up to 0.7
r_num=50 #in a 50 beded ward full of R 
timestep=3
pi_ssr=0.05
abx=0.7
pi_ssr = 1-(1-pi_ssr)^(1/timestep)# 0.03451062 when pi_ssr=0.1, timestep=3
(prob_r = 1-((1-pi_ssr)^r_num)) #0.6513216 when pi_ssr=0.1, timestep=3
(abx = 1-(1-abx)^(1/timestep)) #0.330567 when pi_ssr=0.1, timestep=3

#########mu_r- Haggai Bar-Yoseph, JAC, 2016, 
##colonization 100 -> 76.7% (95% CI=69.3%–82.8%) 1 month ->35.2% (95% CI=28.2%–42.9%) 12 months 
##decolonized 17.2-30.7 1 month, 57.1-71.8 12 months 
t=1:40
mu_r=0.002
#at time t, the probability of s-->S 
becoming_s_at_time_t= (1-mu_r)^(t-1)*mu_r
cumsum_becoming_s_at_time_t=cumsum(becoming_s_at_time_t)
plot(t, cumsum_becoming_s_at_time_t)
#mur=0.002 for ~ 30 days for cummulative risk to become 0.05
#mur=0.02 for ~ 30 days for cummulative risk to become 0.45

####Model 3 
######## K 
mean(d$ini_16S_log) #total capacity in log, K
sd(d$ini_16S_log)

######## r_prop
summary(d$ini_CTXm_copies/d$ini_16S_copies) #proportion of R in total 
#                                            number of Enterbacterobacteriaceae, r_prop
#                                            = exponential distribution 
hist(d$ini_CTXm_copies/d$ini_16S_copies,breaks =seq(0,6,0.01))

####### r_growth, s_growth
#Shaw ISME 2019- microbiome recovery
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
D=1.34
A=-.03
phi1=-1.53
phi2=0.09
y=((D*exp(phi1)*exp(phi2))/(exp(phi2)-exp(phi1)))*(exp(-exp(phi1*t))-exp(-exp(phi2*t)))+A*(1-exp(-exp(phi1*t)))
plot(t,y) #similar distribution os logistic growth 

#####Gibson 2018: doubling time of e coli 15h 
#G (generation time) = (time, in minutes or hours)/n(number of generations), G = t/n
log10(2)/(15/24) #15 hours - a daily growth rate of 115.2 with doubling time 0.625 days
log10(2)/(25/24) #25 hours - a daily growth rate of 0.29 
log10(2)/(5/24)  #15 hours - a daily growth rate of 1.44

t=1:500
r_growth=0.29
K=12
total_capacity=log(0.05*(exp(K)))
r=200
s=1000
r_growth*r*(1 - (r + s)/exp(total_capacity)) #growth in 1 day 
