## summary of the lit review for ABM paper 


d.all = read.csv('metaanalysis/lit_review/abxdur_durtrialappriase_litreview - data.csv')
d = d[which(d$year >= 2000),]
nrow(d)

table(d$interven_type)
12+148+27

59/187

7/187

d = data.frame(x1 = c(7, 18, 4, 18, 8),
              n1 = c(94, 222, 130, 127, 30), 
              x2 = c(19, 38, 12, 17, 5), 
              n2 = c(101, 233, 130, 127, 30))
apply(d, 1, function(x){
  
  out = prop.test(x = c(x['x1'], x['x2']), n = c(x['n1'], x['n2']))
  
  return(c(estimate = round(diff(out$estimate),2)*100, 
         lower.upper = round(out$conf.int,2)*100))
  
})
