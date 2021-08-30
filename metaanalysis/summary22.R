# summary of 19 studies identified 

d = read.csv('metaanalysis/lit_review/clean_extracteddata.csv')
length(unique(d$PMID)) # 24 studies 

only.GP = list()
for (i in unique(d$PMID)){
  dsub = unique(d$bacteria_grp[which(d$PMID == i)])
  if (length(dsub) == 1) {
    if (dsub == 'Gram positives'){
      only.GP[as.character(i)] = i
    }
  }
}

d = d[-which(d$PMID %in% unlist(only.GP)),]
length(unique(d$PMID)) # 22 studies 

selectedPMID = unique(d$PMID)


d = read.csv('metaanalysis/lit_review/extracted_data26.csv')
d = d[which(d$PMID %in% selectedPMID),]

unique(d$country) # only Nigeria low-middle income 
table(d$study_type[!duplicated(d$PMID)])
unique(d$title)
# children: n = 6
# Shortened Antimicrobial Treatment for Acute Otitis Media in Young Children
# Compared to placebo, long-term antibiotics resolve otitis media with effusion (OME) and prevent acute otitis media with perforation (AOMwiP) in a high-risk population: A randomized controlled trial
# Ceftibuten versus trimethoprim-sulfamethoxazole for oral treatment of febrile urinary tract infection in children
# Single-dose extended-release azithromycin versus a 10-day regimen of amoxicillin/clavulanate for the treatment of children with acute otitis media
# Increased risk of acquisition and transmission of ESBL-producing Enterobacteriaceae in malnourished children exposed to amoxicillin
# Meropenem vs standard of care for treatment of neonatal late onset sepsis (NeoMero1): A randomised controlled trial

table(d$setting[!duplicated(d$PMID)])

# intervention arms 
d.int = d[!is.na(d$abx_type),]
d.int$dur.txt = paste(d.int$PMID, d.int$abx_dur)
d.int$abx_dur[!duplicated(d.int$dur.txt)]
summary(as.numeric(d.int$abx_dur[!duplicated(d.int$dur.txt)]))

# change before and after 
bfaft = list()
for (i in selectedPMID){
  
  dsub = d[which(d$PMID == i),]
  pmid = as.character(i)
  
  for (r in unique(dsub$resistance_type)){
    for (b in unique(dsub$bacteria_reported)){
      for (s in unique(dsub$colonisation_site)){
        
        dsubsub = dsub[which(dsub$resistance_type == r & dsub$bacteria_reported == b & dsub$colonisation_site == s),]
        
        if (nrow(dsubsub) > 0) {
          
          abx.long = as.character(max(as.numeric(dsubsub$abx_dur)))
          abx.short = as.character(min(as.numeric(dsubsub$abx_dur)))
          
          pos =  c(dsubsub$n_ind_outcome[which(dsubsub$abx_dur == abx.short)], dsubsub$n_ind_outcome[which(dsubsub$abx_dur == abx.long)])
          denom = c(dsub$n_ind_contributedsamples[which(dsubsub$abx_dur == abx.short)], dsub$n_ind_contributedsamples[which(dsubsub$abx_dur == abx.long)])
          out = prop.test(pos, denom)
          
          bfaft[[pmid]][[paste0(r, b, s)]][['bug']] = paste0(r, b)
          bfaft[[pmid]][[paste0(r, b, s)]][['diff']] = out$estimate[[2]] - out$estimate[[1]] # after - before 
          bfaft[[pmid]][[paste0(r, b, s)]][['estimate']] = out$estimate
          bfaft[[pmid]][[paste0(r, b, s)]][['p']] = out$p.value
        }
      }
    }
  }
}

length(table(unlist(sapply(bfaft, function(x){
  sapply(x, `[[`, "bug")
})))) #types of resistant organisms reported 

length(unlist(sapply(bfaft, function(x){
  sapply(x, `[[`, "diff")
}))) # total number of data points 

summary(unlist(sapply(bfaft, function(x){
  sapply(x, `[[`, "diff")
})))

sum(unlist(sapply(bfaft, function(x){
  sapply(x, `[[`, "p")
})) < 0.05, na.rm = T)

