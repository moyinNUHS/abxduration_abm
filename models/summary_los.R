#frequency summary of patient.matrix - patient id against number of days of stay for each patient

summary_los <- function(patient.matrix){  
  
  # Summarize how often each patient.id is encountered to get days for each id
  los.dur = table(patient.matrix) ##SLOW
  los_duration = array(dim = c(2, max(patient.matrix)))
  # Attach patient ID on 1st row
  los_duration[1,] = 1:max(patient.matrix)
  # Put summary of days on 2nd row
  los_duration[2,] = los.dur
  
  return(los_duration)
}