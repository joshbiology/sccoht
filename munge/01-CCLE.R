#01-CCLE
cache("ccle_rpkm", {
  ccle_rpkm <- readr::read_tsv(file.path(".", "data", "raw", "ccle_rpkm.csv"),col_names = T) %>% 
    df.to.mat
})
