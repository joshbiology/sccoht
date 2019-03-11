#SetupProject.R

### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("ProjectTemplate",
           "doMC",
           "devtools",
           "rfigshare",
           "magrittr",
           "readr",
           "stringr",
           "plyr",
           "tibble",
           "tidyr",
           "dplyr",
           "ggplot2",
           "cowplot",
           "gridExtra",
           "scales",
           "doMC",
           "RColorBrewer",
           "matrixStats",
           "shiny",
           "igraph",
           "ggraph",
           "ggiraph",
           "ggrepel",
           "ggridges",
           "ggsignif",
           "gplots",
           "pheatmap",
           "tibble"
)

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
  install.packages(pkgs = setdiff(pkgs, installed_pkgs))
}


bioc_pkgs <- c("HD2013SGI",
               "diffloop",
               "TxDb.Hsapiens.UCSC.hg19.knownGene",
               "ChIPseeker"
)

#Install bioconductor packages for V 3.8 / R 3.5. This may differ if you're 
#running an older version of R.


if (length(setdiff(bioc_pkgs, installed_pkgs)) > 0) {
  purrr::map(setdiff(pkgs, installed_pkgs), ~ BiocManager::install(., version = "3.8"))
}

### Create directories

dir.create(file.path(".", "output"))
dir.create(file.path(".", "output", "figure1"))
dir.create(file.path(".", "output", "figure2"))
dir.create(file.path(".", "output", "figure3"))
dir.create(file.path(".", "output", "figure4"))
dir.create(file.path(".", "output", "figure5"))
dir.create(file.path(".", "output", "figure6"))
dir.create(file.path(".", "output", "suppfigure1"))
dir.create(file.path(".", "output", "suppfigure2"))
dir.create(file.path(".", "output", "suppfigure3"))
dir.create(file.path(".", "output", "suppfigure5"))

### Download data from figshare
#Only needs to be run once

library(rfigshare)
library(magrittr)
library(purrr)
figshare_meta <- tibble::tribble(
  ~Dataset, ~Figshare_ID,  ~Destination,
  "remaining_data_id", 7801718, ".")

download_figshare_record <- function(record_id, data_dir) {
  figshare_article <- fs_details(record_id)
  
  for (figshare_file in figshare_article$files) {
    cat("downloading", figshare_file$name, "\n")
    file_path <- file.path(data_dir, figshare_file$name)
    download.file(figshare_file$download_url, file_path)
  }
}

walk2(figshare_meta$Figshare_ID,
     figshare_meta$Destination, function(x,y) download_figshare_record(x,y))

library(ProjectTemplate)
load.project(munging = T)

sessionInfo()

#Run on Josh's local computer after 3-11-2019 data / code freeze
# R version 3.5.1 (2018-07-02)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS High Sierra 10.13.6
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats4    parallel  grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] bindrcpp_0.2.2        pheatmap_1.0.10       lubridate_1.7.4       reshape2_1.4.3        dependr_1.1          
# [6] taigr_2.6.1.8         httr_1.3.1            RCurl_1.95-4.11       bitops_1.0-6          RColorBrewer_1.1-2   
# [11] purrr_0.2.5           scales_1.0.0.9000     gridExtra_2.3         cowplot_0.9.3         ggplot2_3.1.0        
# [16] dplyr_0.7.8           tidyr_0.8.2           tibble_1.4.2          plyr_1.8.4            ChIPpeakAnno_3.16.1  
# [21] VennDiagram_1.6.20    futile.logger_1.4.3   GenomicRanges_1.34.0  GenomeInfoDb_1.18.1   Biostrings_2.50.2    
# [26] XVector_0.22.0        IRanges_2.16.0        S4Vectors_0.20.1      BiocGenerics_0.28.0   stringr_1.3.1        
# [31] readr_1.2.1           magrittr_1.5          usethis_1.4.0         devtools_2.0.1        ProjectTemplate_0.8.2
# 
# loaded via a namespace (and not attached):
#   [1] ChIPseeker_1.18.0                       backports_1.1.2                         fastmatch_1.1-0                        
# [4] igraph_1.2.2                            lazyeval_0.2.1                          splines_3.5.1                          
# [7] BiocParallel_1.16.5                     gridBase_0.4-7                          urltools_1.7.2                         
# [10] digest_0.6.18                           foreach_1.4.4                           ensembldb_2.6.3                        
# [13] GOSemSim_2.8.0                          viridis_0.5.1                           GO.db_3.7.0                            
# [16] gdata_2.18.0                            memoise_1.1.0                           BSgenome_1.50.0                        
# [19] limma_3.38.2                            remotes_2.0.2                           matrixStats_0.54.0                     
# [22] enrichplot_1.2.0                        prettyunits_1.0.2                       colorspace_1.3-2                       
# [25] blob_1.1.1                              ggrepel_0.8.0                           jsonlite_1.5                           
# [28] callr_3.0.0                             crayon_1.3.4                            Sushi_1.20.0                           
# [31] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 graph_1.60.0                            bindr_0.1.1                            
# [34] survival_2.43-1                         zoo_1.8-4                               iterators_1.0.10                       
# [37] glue_1.3.0                              gtable_0.2.0                            zlibbioc_1.28.0                        
# [40] UpSetR_1.3.3                            seqinr_3.4-5                            DelayedArray_0.8.0                     
# [43] pkgbuild_1.0.2                          DOSE_3.8.2                              futile.options_1.0.1                   
# [46] DBI_1.0.0                               edgeR_3.24.3                            Rcpp_1.0.0                             
# [49] plotrix_3.7-4                           viridisLite_0.3.0                       progress_1.2.0                         
# [52] units_0.6-2                             annotables_0.1.91                       gridGraphics_0.3-0                     
# [55] europepmc_0.3                           bit_1.1-14                              fgsea_1.8.0                            
# [58] gplots_3.0.1                            pkgconfig_2.0.2                         XML_3.98-1.16                          
# [61] farver_1.1.0                            locfit_1.5-9.1                          ggplotify_0.0.3                        
# [64] tidyselect_0.2.5                        rlang_0.3.0.1                           AnnotationDbi_1.44.0                   
# [67] munsell_0.5.0                           tools_3.5.1                             cli_1.0.1                              
# [70] moments_0.14                            RSQLite_2.1.1                           ade4_1.7-13                            
# [73] ggridges_0.5.1                          yaml_2.2.0                              processx_3.2.0                         
# [76] bit64_0.9-7                             fs_1.2.6                                caTools_1.17.1.1                       
# [79] AnnotationFilter_1.6.0                  ggraph_1.0.2                            pbapply_1.3-4                          
# [82] RBGL_1.58.1                             formatR_1.5                             xml2_1.2.0                             
# [85] DO.db_2.9                               biomaRt_2.38.0                          compiler_3.5.1                         
# [88] rstudioapi_0.8                          curl_3.2                                statmod_1.4.30                         
# [91] tweenr_1.0.1                            stringi_1.2.4                           idr_1.2                                
# [94] ps_1.2.1                                GenomicFeatures_1.34.1                  desc_1.2.0                             
# [97] lattice_0.20-38                         ProtGenerics_1.14.0                     Matrix_1.2-15                          
# [100] multtest_2.38.0                         pillar_1.3.0                            BiocManager_1.30.4                     
# [103] triebeard_0.3.0                         data.table_1.11.8                       qvalue_2.14.1                          
# [106] rtracklayer_1.42.1                      R6_2.3.0                                KernSmooth_2.23-15                     
# [109] sessioninfo_1.1.1                       codetools_0.2-15                        lambda.r_1.2.3                         
# [112] gtools_3.8.1                            boot_1.3-20                             MASS_7.3-51.1                          
# [115] assertthat_0.2.0                        pkgload_1.0.2                           SummarizedExperiment_1.12.0            
# [118] rprojroot_1.3-2                         withr_2.1.2                             regioneR_1.14.0                        
# [121] diffloop_1.10.0                         GenomicAlignments_1.18.1                Rsamtools_1.34.0                       
# [124] GenomeInfoDbData_1.2.0                  hms_0.4.2                               rvcheck_0.1.3                          
# [127] ggforce_0.1.3                           Biobase_2.42.0                          base64enc_0.1-3    