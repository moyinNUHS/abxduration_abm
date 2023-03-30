library(metafor)

### Enter data from 5 trials 
dat = matrix(c(1, 7, 8, 272, 7, 7, 94, 19, 19, 101, 11, 22, 1, 
               2, 5, 10, 520, 13, 28, 222, 16, 38, 233, 4, 11, -3,
               3,	1,	5,	260,	9,	7,	77,	18,	12,	65,	9,	22,	-3,
               4,	3,	5,	324,	12,	8,	66,	19,	12,	63,	7,	21,	-7,	
               5,	3,	14,	60,	27,	8,	30,	17,	5,	30,	-10,	4,	-14
), nrow = 5, byrow = T)
dat = as.data.frame(dat)
colnames(dat) = c('trial_no', 'short_dur', 'long_dur', 'n_participant', 'resprop_short', 'resn_short', 'n_short', 
                  'resprop_long', 'resn_long', 'n_long', 'diff', 'diff_upp', 'diff_low')
dat$author = c('Lutsar', 'Hoberman', 'Ceran', 'Merode', 'Dow')
dat$year = c(2020, 2016, 2010, 2005, 2004)
dat$short_pos = dat$resn_short
dat$long_pos = dat$resn_long
dat$short_neg = dat$n_short - dat$resn_short
dat$long_neg = dat$n_long - dat$resn_long

### calculate log risk ratios and corresponding sampling variances
out <- escalc(measure="RD", ai = long_pos, bi = long_neg, ci = short_pos, di = short_neg, data = dat)

### fit random-effects model
res <- rma(yi, vi, data = out, slab = paste(author, year, sep=", "))

png(file = "~/Documents/nBox/angelsfly/indiv_abxduration/manuscript/manuscript/graphs/funnel.png", 
    width = 1000, # The width of the plot in inches
    height = 600) # The height of the plot in inches
funnel(res, yaxis="sqrtni", label="out", xlim = c(-0.2, 0.2))
dev.off()
