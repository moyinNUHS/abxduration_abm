## summary of the lit review for ABM paper 


d.all = read.csv('metaanalysis/lit_review/abxdur_durtrialappriase_litreview - data.csv')
d = d[which(d$year >= 2000),]
nrow(d)

table(d$interven_type)
12+148+27

59/187

7/187

d = data.frame(x1 = c(5, 28, 7, 8, 8), # short outcome
               n1 = c(70, 222, 77, 64, 30), #short n
               x2 = c(5, 38, 12, 9, 5), #long outcome
               n2 = c(124, 233, 65, 60, 30)) # long n
apply(d, 1, function(x){
  
  out = prop.test(x = c(x['x1'], x['x2']), n = c(x['n1'], x['n2']))
  
  return(c(estimate = round(diff(out$estimate),2)*100, 
           lower.upper = sort(-round(out$conf.int,2)*100)))
  
})
