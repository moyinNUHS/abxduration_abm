##########################################################################
#######Effect of antibiotic duration on hospitalised patients ############
############## Show individual vs population dynamics   ##################
##########################################################################

# clean working environment
rm(list=ls()) 

# load libraries 
library(ggplot2)
library(ggpubr)
library(reshape)

source('models/get_output_indivVSpop_populationgrowth.R')

# K = exp(20)
# total_prop = 0.1
# p.infect = 0.4
prop_R = 0.25
# pi_ssr = 0.5
# r_trans = 0.5
# fitness.r = 1
# r_thres = 0.05
# s_growth = 1
# abx.s = 0.7
# abx.r = 0.7
# max.los = 14
max.dur = 21

out = run_indivVSpop_populationgrowth(K = exp(20),
                                      max.los = 14,
                                      #no antibiotic prescribed after admission
                                      p.infect.after = 0.00001,
                                      p.r.day1 = 1,
                                      p.r.after = 0.5, #does not matter
                                      total_prop = 0.5,
                                      p.infect = 0.5,
                                      prop_R =  prop_R,
                                      pi_ssr = 0.7,
                                      r_trans = 0.7,
                                      fitness.r = 1,
                                      r_thres = 0.01,
                                      s_growth = 0.5,
                                      abx.s = 0.5,
                                      abx.r = 0.5,
                                      iterations = 20,
                                      max.dur = max.dur)

#saveRDS(out, file = 'runs/indivVSpop.RDS')

#out = readRDS('runs/indivVSpop.RDS')

## within host 
eff.lab = 'Administered antibiotics\nactive against susceptible\nand resistant organisms'
ineff.lab = 'Administered antibiotics\nactive only against\nsusceptible organisms'
dw = out$withinhost
dw$abx<- factor(dw$abx, levels = c("eff", "ineff"),
                labels = c(eff.lab,
                           ineff.lab))

w = ggplot(dw, aes(x = day, y = value, color = res)) + 
  geom_line(size = 1.2) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(y = "Number of bacteria \n(proportion of max gut capacity)", 
       x = 'Day of antibiotic treatment',
       title = 'Effect of antibiotic on the treated individuals') +
  scale_color_manual(values = c('#A31621', '#12664F'), 
                     name = ' \n ',
                     labels = c('Resistant bacteria','Susceptible bacteria')) + 
  facet_grid(abx~.) + 
  theme_minimal() + 
  theme(legend.position = 'bottom', 
        text = element_text(size = 17),
        plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        strip.text.y = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


##population 
dp = out$pop
dp$abx<- factor(dp$abx, levels = c("Effective antibiotic available", "Effective antibiotic not available"),
                labels = c(eff.lab,
                           ineff.lab))
dp$variable = factor(dp$variable, levels = rev(levels(dp$variable)))

line.dat = data.frame(y = c( prop_R,  prop_R), x = c(0.5, max.dur-1), variable = c('Transmitted'))
r.ineff = sum(dp$value[which(dp$abx == ineff.lab & dp$Day == (max.dur-1) & dp$variable != 'Non-carrier')]/50)
r.eff = sum(dp$value[which(dp$abx ==  eff.lab & dp$Day == (max.dur-1) & dp$variable != 'Non-carrier')]/50)
arrow.dat = data.frame(y = c(0, r.eff , r.eff , 1, 0, r.ineff, 1, r.ineff ),
                       yend = c(r.eff , 0, 1, r.eff , r.ineff, 0, r.ineff, 1),
                       x = rep(max.dur, 8),
                       xend = rep(max.dur, 8), 
                       variable = c('Transmitted'), 
                       abx = c(rep(ineff.lab, 4), 
                               rep(eff.lab, 4)))
lab.dat = data.frame(y = c(0.6, 0.1, 0.65, 0.155),
                     x = c(max.dur+1.5, max.dur+1.5, max.dur+1.5, max.dur+1.5),
                     variable = c('Transmitted'), 
                     abx = c(rep(ineff.lab, 2), 
                             rep(eff.lab, 2)),
                     lab = c('Non-carriers', 'Resistance\ncarriers', 
                             'Non-carriers', 'Resistance\ncarriers'))


p = ggplot(dp, aes(x = Day, fill = variable, y = value)) +
  geom_bar(position="fill", stat="identity", alpha = 0.7, width = 1) + 
  labs(y = 'Proportion of patients in the ward', 
       x = 'Mean duration of antibiotic treatment received', 
       title = 'Effect of antibiotic on the ward population') +
  geom_text(aes(x = 5, y =  prop_R + 0.085, label = 'Prevalence of resistance\ncarriers at admission'), color = 'grey10', size = 4.5) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  #scale_x_continuous(breaks = seq(1, 20, by = 5)) +
  scale_fill_manual(values = rev(c('#D88C9A', '#931621', '#f7e6e7', '#BAD9B5')), 
                    labels = c('Selection', 'Transmission', 'Admitted as R carrier'),
                    breaks = c('Selection', 'Transmitted', 'Admitted as R'),
                    name = c('Type of resistance\nacquisition')) +
  facet_grid(abx~.) + 
  theme_minimal() + 
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        text = element_text(size = 17), 
        axis.text.y=element_blank(),
        axis.line.x = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.text.y = element_text(size = 17, face = 'bold')) + 
  geom_line(aes(x=x, y=y), dat = line.dat, linetype = 'dashed', color = 'grey10', size = 1.2) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), arrow = arrow(length = unit(0.25, "cm")), data = arrow.dat) + 
  geom_text(aes(x = x, y = y, label = lab), dat = lab.dat, color = 'Black', size = 5, angle = 270) + 
  guides(fill = guide_legend(nrow = 2))
  

ggarrange(w, p, ncol = 2)  
ggsave('~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/final_main/indivVSpop.png',
       width = 14,
       height = 10)

