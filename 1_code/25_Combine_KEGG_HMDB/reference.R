match_compound_databases <- 
  function(database1, 
           database2, 
           match_id = "KEGG_ID"){
  id1 <- database1 %>% 
    dplyr::pull(match_id)
  id2 <-
    database2 %>% 
    dplyr::pull(match_id)
  
  match_idx <-
    match(id2, id1)
  
  database1_new <-
  seq_len(match_idx) %>% 
  purrr::map(function(i){
    cat(i, " ")
    if(!is.na(match_idx[i])){
      ###Synonyms1
      Synonyms1 <- 
        database1[match_idx[i],]$Synonyms %>% 
        stringr::str_split("\\{\\}") %>%
        `[[`(1) %>% 
        unique()
        
        
        Synonyms2 <- 
        c(database2[i,]$Compound_name,
        stringr::str_split(database2[i,]$Synonyms, "\\{\\}")[[1]]) %>%  
        unique()
        
        Synonyms2 <-
          Synonyms2[!is.na(Synonyms2)]
        
        Synonyms <-
        paste0(unique(c(Synonyms1, Synonyms2)), collapse = "{}")
        
        database1[match_idx[i],]$Synonyms <- Synonyms
        
        ##
        
    }
  })
  
  
  database2_new <-
    database2[which(is.na(match_idx)),]

  
}