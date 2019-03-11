### 06-Differential_Expression
clean_deseq_output <- function(path, name) {
	read_csv(path) %>% 
		dplyr::rename(Gene = X1) %>% 
		filter(!is.na(padj)) %>% 
		dplyr::select(-lfcSE, -stat) %>% 
		mutate(Condition = name)
}
	
# Differential expression
#SCCOHT
path <- file.path(".", "data", "raw", "BIN67_DESeq")
filenames <- list.files(path, ".csv")
bin_deseq <- adply(filenames, 1, 
									 function(x) clean_deseq_output(file.path(path, x), 
									 															 x %>% str_sub(end = -11))) %>%
	as.tibble() %>% 
	dplyr::select(-1)


#DESEQ complete
cache("bin_deseq_complete", {
	#Make version of bin_deseq in which fold changes are calculated for all genes
	bin_gexp_mat <- bin_gexp %>% 
		select(-Cell_Line) %>% 
		filter(Condition != "SMARCA2") %>% 
		group_by(Gene, Condition) %>% 
		summarize(Mean_RPKM = mean(RPKM)) %>% 
		df.to.mat()
	
	bin_gexp_mat <- bin_gexp_mat + 1 #RPKM + 1
	
	bin_gexp_fc <- log(bin_gexp_mat/bin_gexp_mat[,"Control"], 2) %>% 
		mat.to.df(row.name = "Gene", col.name = "Condition", dat.name = "Log2FC")
	
	#Make version of DESEQ in which missing values are added as 0
	tmp_bin_deseq_mat <- bin_deseq %>% 
		select(Gene, Condition, log2FoldChange) %>% 
		df.to.mat()
	
	tmp_bin_deseq_mat[is.na(tmp_bin_deseq_mat)] <- 0
	
	bin_deseq_complete <- tmp_bin_deseq_mat %>% 
		mat.to.df(row.name = "Gene", col.name = "Condition", dat.name = "log2FoldChange")
})

#Up and Down Genesets
bin_wt_up <- bin_deseq %>% 
	filter(Condition == "SCCOHT_SMARCA4") %>% 
	filter(log2FoldChange > 2, 
				 padj < 0.05,
				 Gene %in% bin_expressed) %>% 
	pull(Gene)

bin_wt_down <- bin_deseq %>% 
	filter(Condition == "SCCOHT_SMARCA4") %>% 
	filter(log2FoldChange < -2, 
				 padj < 0.05,
				 Gene %in% bin_expressed) %>% 
	pull(Gene)
