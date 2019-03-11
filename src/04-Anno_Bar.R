#04-Anno_Bar.R

library(ProjectTemplate)
load.project()
source(file.path(".", "munge", "05-Gene_Expression.R"))
source(file.path(".", "munge", "07-MACS_peaks.R"))

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#Figure 2G
c(single_epitopes[3:4], complex_k785r[["Ctrl_Residual_Reference_Peaks"]]) %>% 
  set_names(c("PBAF_SMARCA4", "BAF_SMARCA4", "Z_Combined_residual")) %>% 
  lapply(ChIPseeker::annotatePeak, TxDb=txdb,
         tssRegion=c(-3000, 3000), verbose=FALSE) %>% 
  plotAnnoBar() +
  ggsave(file.path(".", "output", "figure2", "fig2g_anno_bar.pdf"), width = 6, height = 3)