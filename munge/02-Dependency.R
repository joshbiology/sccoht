#02-Dependency
cache("avana_dep", {
	gene_effect <- read_csv(file.path(".", "data", "raw", "gene_effect.csv"))
	avana_dep <- gene_effect[-1]  %>%
		as.matrix()
	
	rownames(avana_dep) <- gene_effect[[1]]
	colnames(avana_dep) <- colnames(avana_dep) %>% word()
	avana_dep <- avana_dep %>% t()
})

cache("avana_analysis_genes", {
  source(file.path(".", "munge", "01-CCLE.R"))
  globalSD <- sd(avana_dep, na.rm=T)
  most_dependent_line <- adply(avana_dep, 1, function(x) {
    xmin <- which.min(x)
    data_frame(CellLine = names(xmin),
               Dep = x[xmin],
               SD = ((x[xmin] - mean(x, na.rm=T)) / globalSD))
  }, .id="Gene")
  
  xpr <- most_dependent_line %>%
    select(CellLine, Gene) %>%
    filter(CellLine %in% colnames(ccle_rpkm),
           Gene %in% rownames(ccle_rpkm)) %>%
    mutate(RPKM = ccle_rpkm[as.matrix(.[,c("Gene", "CellLine")])])
  
  analysis_genes <- most_dependent_line %>%
    left_join(xpr) %>%
    filter(Dep < -0.3 & (is.na(RPKM) | RPKM > 0))
  return(as.character(analysis_genes$Gene))
})