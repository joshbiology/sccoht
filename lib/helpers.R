#helpers.R

bed_subtract <- function(A, B) {
  A[!A %over% B]
}


#Peak to Gene mapping
assign_peakset_to_geneset <- function(lst, tss_set) {
  lst %>% 
    llply(function(x) {
      findOverlaps(tss_set, x) %>% 
        as_tibble() %>% 
        dplyr::count(queryHits)
    }) %>% 
    ldply() %>% 
    mutate(Gene = tss_set$Gene[queryHits]) %>% 
    filter(Gene %in% bin_expressed) %>% 
    select(Gene, Peakset = .id, Count = n) %>% 
    unique() %>% 
    group_by(Gene) %>% 
    mutate(Total = sum(Count)) %>% 
    arrange(-Count) %>% 
    filter(row_number()==1) %>% 
    ungroup()
}

assign_peakset_to_ref_set <- function(lst, ref_set) {
  lst %>% 
    llply(function(x) {
      findOverlaps(ref_set, x) %>% 
        as_tibble() %>% 
        dplyr::count(queryHits)
    }) %>% 
    ldply() %>% 
    as.tibble() %>% 
    select(queryHits, Peakset = .id, Count = n) %>% 
    unique() %>% 
    group_by(queryHits) %>% 
    mutate(Total = sum(Count)) %>% 
    arrange(-Count) %>% 
    filter(row_number()==1) %>% 
    ungroup() 
}

custom_fishers_greater <- function(a,b,c,d) {
  #https://stackoverflow.com/questions/29214222/error-all-entries-of-x-must-be-nonnegative-and-finite-in-fisher-test
  fisher.test(matrix(c(a, b, c, d) %>% unlist(), nrow = 2),
              alternative = "greater")
}