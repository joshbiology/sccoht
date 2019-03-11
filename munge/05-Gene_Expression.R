### 05-Gene_Expression

# BIN-67 lentiviral reintroductions
bin_gexp <- file.path(".", "data", "raw", "BIN67_RNAseq") %>% 
	file.path(., list.files(., pattern = ".txt")) %>% 
	lapply(read_tsv) %>% 
	Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2), .) %>% 
	gather(Cell_Line, RPKM, -1) %>% 
	separate(Cell_Line, into = c("Cell_Line", "Condition", "Replicate"), sep = "_") %>% 
	group_by(Cell_Line, Condition, Gene) %>% 
	mutate(Mean_RPKM = mean(RPKM)) %>% 
	ungroup()

#Total expressed genes (used in downstream analysis)
expression_threshold <- 2 #RPKM

bin_expressed <- bin_gexp %>% 
  filter(RPKM > expression_threshold) %>% 
  pull(Gene) %>% 
  unique() #13189

# Tumor data from Nat Gen 2015 SMARCA4-DTS paper (code credit to John Pulice)
cache("tumor_gexp", {
	path <- "./data/raw/le_loarer_2015_PMID_ 26343384/"
	countfiles <- list.files(path, pattern=".*read_cnt$")
	countnames <- sapply(strsplit(countfiles, split="[.]"), function(x) x[1])
	countfiles <- countfiles[order(countnames)]
	countnames <- countnames[order(countnames)]
	
	meta <- cbind(countnames, sapply(strsplit(countnames, split="_"), function(x) x[1]),sapply(strsplit(countnames, split="_"), function(x) x[2]))
	
	for(i in seq_len(length(countfiles))){
		temp <- read.table(file.path(path, countfiles[i]))
		if(i == 1){
			counts <- array(0,c(nrow(temp),length(countfiles)))
			rownames(counts) <- temp[,1]
			colnames(counts) <- countnames
		}
		counts[,i] <- temp[,5]
	}
	
	tumor_gexp <- counts[, meta[,3] == "Tumor" & meta[,2] != "LCa"] %>% 
		as.tibble() %>% 
		add_column(Gene = rownames(counts), .before = 1) %>% 
		gather(Class, RPKM, -1) %>% 
		separate(Class, into = c("Class", "Type", "Sample"), sep = "_")
})




