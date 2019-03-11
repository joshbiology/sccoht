### 07-MACS peaks

#TODO: Rename luciferase to Control in the relevant file names
#TODO: Delete unused peak files

clean_macs2_output <- function(path, peak_type = "narrowPeak", peak_names = NA) {
	# Import all peak files in that folder
	peak_files <- list.files(path,
													 pattern=paste(".*", peak_type, "$", sep = ""))
	
	# Optional: derive list names from file names
	if (is.na(peak_names)) {
		peak_names <- sapply(strsplit(peak_files, split="_peaks"), function(x) x[1])
	}
	
	# Create list of peaks (GRangesList)
	peak_files %>% 
		alply(1, function(x) diffloop::bedToGRanges(file.path(path, x))) %>% 
		set_names(peak_names) %>% 
		GRangesList() %>% 
		llply(granges) %>% 
		llply(fix_map)
}

#Fixing and mapping peaks
chrom <- paste0("chr",c(1:22,"X","Y"))
blacklist <- diffloop::bedToGRanges(file.path(".", "data", "raw", "wgEncodeDacMapabilityConsensusExcludable.bed"))

fix_map <- function(temp) {
	temp <- temp %>% 
		bed_subtract(blacklist)
	temp[seqnames(temp) %in% chrom]
}

#SCCOHT
cache("bin_chip_peaks", {
	bin_chip_peaks <- clean_macs2_output(file.path(".", "data", "raw", "chip_peaks"))
	
})

cache("bin_broad_peaks", {
	bin_broad_peaks <- clean_macs2_output(file.path(".", "data", "raw", "chip_peaks"),
																				peak_type = "broadPeak")
	
})

# TSS "peaks"

generate_tss <- function(flank) {
	annotables::grch37 %>% 
		mutate(TSS = ifelse(strand == 1, start, end)) %>% 
		mutate(start_promoter = ifelse(strand == 1, start - flank/2, end - flank/2),
					 end_promoter = ifelse(strand == 1, start + flank/2, end + flank/2)) %>% 
		mutate(strand = case_when(strand ==  1 ~ "+",
															strand == -1 ~ "-")) %>% 
		select(Gene = symbol, TSS, start = start_promoter, end = end_promoter, chr, strand) %>%
		filter(Gene %in% bin_expressed) %>% 
		makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
		diffloop::addchr() %>% 
		fix_map()
}

tss_2kb <- generate_tss(2000)
tss_4kb <- generate_tss(4000)
tss_5kb <- generate_tss(5000)
tss_30kb <- generate_tss(30000)



# Reference Peaks

read_peaks_from_path <- function(path) {
  purrr::map(list.files(path, full.names = T),
             ChIPseeker::readPeakFile) %>% 
    purrr::set_names(list.files(path, full.names = F) %>% word(1, sep = "\\."))
}

bin_atac_peaks <- read_peaks_from_path(file.path(".", "data", "raw", "atac_q1e-3"))

single_epitopes <- read_peaks_from_path(file.path(".", "data", "interim", "single_epitopes"))
complex_k785r <- read_peaks_from_path(file.path(".", "data", "interim", "complex_k785r"))

enh <- IRanges::subsetByOverlaps(bin_broad_peaks[["BIN67_SMARCA4_H3K4me1"]],
                        bin_atac_peaks[["ATAC_SMARCA4_Merged_peaks"]])

