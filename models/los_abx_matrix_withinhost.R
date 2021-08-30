# here we simulate one individual 

# generate a table of number of days we want to observe (rows) -
# against number of repeats for each individual per set of parameters (denoted as n.bed as columns), filled in with repeat id numbers
# stays for the entire observation time 

# differences with los_abx_matrix_varydur: 
# 1. los is the same as the observation time - no max.los
# 2. get antibiotics on day 1, all broad spectrum 


los_abx_table_withinhost <- function(n.bed, n.day, short_dur, long_dur, timestep = 1, sd = 1){
  
  
  patient.matrix = matrix(rep(1:n.bed, each = n.day), ncol = n.bed)
  abx.matrix.short = matrix(rep(rep(c(2,0), c(short_dur, n.day - short_dur)), each = n.bed), ncol = n.bed, byrow = T)
  abx.matrix.long = matrix(rep(rep(c(2,0), c(long_dur, n.day - long_dur)), each = n.bed), ncol = n.bed, byrow = T)
  

  return(list(patient.matrix, abx.matrix.short, abx.matrix.long))
}




