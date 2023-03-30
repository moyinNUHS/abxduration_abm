# screen for duplicates 

setwd('~/Documents/nBox/git_projects/abxduration_abm/metaanalysis/')

get.deduplicated <- function(purpose){
  
  filenames = paste0(paste0('lit_review/raw_files/', purpose, '/'), list.files(paste0('lit_review/raw_files/', purpose)))
  filenames = filenames[grep('oct22', filenames)]
  pubmed.filenames = filenames[grep('pubmed', filenames)]
  embase.filenames = filenames[grep('embase', filenames)]
  
  raw.database.list = lapply(list(pubmed.filenames, embase.filenames), function(x){
    
    csv.list = lapply(x, read.csv)
    combine = do.call('rbind', csv.list)
    
    colnames(combine)[grep('PMID', colnames(combine))] = 'PMID'
    
    if (!'PUI' %in% colnames(combine)) {combine$PUI = NA}
    
    combine[, c('PMID', 'PUI', 'Title')]
  })
  
  raw = do.call('rbind', raw.database.list)
  
  message(paste('There are', nrow(raw.database.list[[1]]), 'total MEDLINE titles.'))
  message(paste('There are', nrow(raw.database.list[[2]]), 'total EMBASE titles.'))
  
  # number of duplicates based on PMID 
  sum(duplicated(raw$PMID))
  raw.unique = raw[!duplicated(raw$PMID),]

  raw.unique = raw.unique[!duplicated(raw.unique$Title),]
  
  message(paste('There are', nrow(raw) , 'total titles.'))
  message(paste('There are', nrow(raw) - nrow(raw.unique), 'duplicates removed from PMID and titles.'))
  
  write.csv(raw.unique, file = paste0('lit_review/unique_entries_', purpose, Sys.Date(), '.csv'))
}

## for duration and colonising bacteria metaanalysis 2000 - 2021 
get.deduplicated(purpose = 'abxdur')

## 1920 - 1999 
get.deduplicated(purpose = 'durtrialappraise')



