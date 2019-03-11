#01_Readcount_Scatters.R

library(magrittr)
library(tidyverse)

##All readcounts cached for ease of plot generation. 
#Important note: all processing done with Rsubread 1.24.2.

#### ChIP-seq scatters ####

c1_reads <- read_tsv(file.path(".", "data", "interim", "rsubread", "c1_reads.tsv")) %>%
  as.matrix()

cor(c1_reads[,1], c1_reads[,2], method = "spear")
cor(c1_reads[,1], c1_reads[,3], method = "spear")
cor(c1_reads[,2], c1_reads[,3], method = "spear")

pdf(file.path(".", "output", "suppfigure1", "sfig1f_control_scatter_1.pdf"), width = 5, height = 5)
smoothScatter(log(c1_reads[,1], 2), 
              log(c1_reads[,2], 2), 
              colramp = colorRampPalette(RColorBrewer::brewer.pal(9, "Greys")),
              xlab = "SMARCC1 rep 1 log2RPKM", 
              ylab = "SMARCC1 rep 2 log2RPKM")
dev.off()

pdf(file.path(".", "output", "suppfigure1", "sfig1f_control_scatter_2.pdf"), width = 5, height = 5)
smoothScatter(log(c1_reads[,1], 2), 
              log(c1_reads[,3], 2), 
              colramp = colorRampPalette(RColorBrewer::brewer.pal(9, "Greys")),
              xlab = "SMARCC1 rep 1 log2RPKM", 
              ylab = "SMARCC1 rep 3 log2RPKM")
dev.off()

pdf(file.path(".", "output", "suppfigure1", "sfig1f_control_scatter_3.pdf"), width = 5, height = 5)
smoothScatter(log(c1_reads[,2], 2), 
              log(c1_reads[,3], 2), 
              colramp = colorRampPalette(RColorBrewer::brewer.pal(9, "Greys")),
              xlab = "SMARCC1 rep 1 log2RPKM", 
              ylab = "SMARCC1 rep 3 log2RPKM")
dev.off()


#Figure 2C - control ChIP-seq
fig2c <- readr::read_tsv(file.path(".", "data", "interim", "rsubread", "control_dpf2_arid2_counts.tsv"))
pdf(file.path(".", "output", "figure2", "fig2c_control_scatter.pdf"), width = 5, height = 5)
smoothScatter(fig2c %>% 
                select(BIN67_Control_200, BIN67_Control_45D) %>% 
                as.matrix() %>% 
                "+"(1) %>% 
                log(2))
dev.off()
cor(fig2c$BIN67_Control_200, fig2c$BIN67_Control_45D, method = "spearman")

#Figure 2E - rescue SMARCA4 ChIP-seq
fig2e <- readr::read_tsv(file.path(".", "data", "interim", "rsubread", "a4_ss18_arid2_counts.tsv"))
pdf(file.path(".", "output", "figure2", "fig2e_smarca4_rescue_scatter.pdf"), width = 5, height = 5)
smoothScatter(fig2e %>% 
                select(BIN67_A4_200, BIN67_A4_SS18) %>% 
                as.matrix() %>% 
                "+"(1) %>% 
                log(2))
dev.off()

#Supplementary Figure 3F - SS18 vs ARID2 rescue T910M/K785R SMARCA4 ChIP-seq
suppfig3f <- readr::read_tsv(file.path(".", "data", "interim", "rsubread", "t910m_ss18_arid2_counts.tsv"))
pdf(file.path(".", "output", "suppfigure3", "suppfig3f_t910m_rescue_scatter.pdf"), width = 5, height = 5)
smoothScatter(suppfig3f %>% 
                select(BIN67_T910M_ARID2, BIN67_T910M_SS18) %>% 
                as.matrix() %>% 
                "+"(1) %>% 
                log(2),
              xlim = c(0,8),
              ylim = c(0,8))
dev.off()

suppfig3fb <- readr::read_tsv(file.path(".", "data", "interim", "rsubread", "k785r_ss18_arid2_counts.tsv"))
pdf(file.path(".", "output", "suppfigure3", "suppfig3f_k785r_rescue_scatter.pdf"), width = 5, height = 5)
smoothScatter(suppfig3fb %>% 
                select(BIN67_K785R_ARID2_ZM101, BIN67_K785R_SS18_ZM103) %>% 
                as.matrix() %>% 
                "+"(1) %>% 
                log(2),
              xlab = "BIN67_K785R_SS18",
              ylab = "BIN67_K785R_ARID2",
              xlim = c(0,8),
              ylim = c(0,8))
dev.off()


#Supplementary Figure 3E - SMARCA4 vs SMARCC1 ChIP-seq over various reintroduction conditions
suppfig3e <- readr::read_tsv(file.path(".", "data", "interim", "rsubread", "rep_a4_c1_counts.tsv"))
pdf(file.path(".", "output", "suppfigure3", "suppfig3e_k785r_rescue_scatter.pdf"), width = 5, height = 5)
smoothScatter(suppfig3e %>% 
                select(K785R_BRG1, K785R_BAF155) %>% 
                as.matrix() %>% 
                "+"(1) %>% 
                log(2),
              xlab = "K785R_BRG1",
              ylab = "K785R_BAF155")
dev.off()

cor(suppfig3e$K785R_BAF155, suppfig3e$K785R_BRG1, method = "spearman")

pdf(file.path(".", "output", "suppfigure3", "suppfig3e_t910m_rescue_scatter.pdf"), width = 5, height = 5)
smoothScatter(suppfig3e %>% 
                select(T910M_BRG1, T910M_BAF155) %>% 
                as.matrix() %>% 
                "+"(1) %>% 
                log(2),
              xlab = "T910M_BRG1",
              ylab = "T910M_BAF155")
dev.off()

cor(suppfig3e$T910M_BRG1, suppfig3e$T910M_BAF155, method = "spearman")

pdf(file.path(".", "output", "suppfigure3", "suppfig3e_wt_rescue_scatter.pdf"), width = 5, height = 5)
smoothScatter(suppfig3e %>% 
                select(SMARCA4_BRG1, SMARCA4_BAF155) %>% 
                as.matrix() %>% 
                "+"(1) %>% 
                log(2),
              xlab = "SMARCA4_BRG1",
              ylab = "SMARCA4_BAF155")
dev.off()

cor(suppfig3e$SMARCA4_BRG1, suppfig3e$SMARCA4_BAF155)



#### ATAC-seq scatters ####

atac_df <- read_tsv(file.path(".", "data", "interim", "rsubread", "atac_counts.tsv"))

atac_df <- atac_df  %>% 
  mutate(T910_FC = (T910M+1)/(Ctrl+1),
         SMARCA4_FC = (SMARCA4+1)/(Ctrl+1),
         SMARCA2_FC = (SMARCA2+1)/(Ctrl+1),
         K785R_FC = (K785R+1)/(Ctrl+1)) 

colnames(atac_df)
stats::cor(atac_df$SMARCA4_FC, atac_df$SMARCA2_FC, use = "pairwise.complete.obs")
stats::cor(atac_df$SMARCA4_FC, atac_df$T910_FC, use = "pairwise.complete.obs")
stats::cor(atac_df$SMARCA4_FC, atac_df$K785R_FC, use = "pairwise.complete.obs")

lims <- c(0.5, 3) 
pdf(file.path(".", "output", "figure3", "fig3b_wt_rescue_scatter_wt.pdf"), width = 5, height = 5)
smoothScatter(log(atac_df$SMARCA4_FC+1, 2), 
              log(atac_df$SMARCA2_FC+1, 2), 
              colramp = colorRampPalette(RColorBrewer::brewer.pal(9, "Greys")),
              xlim = lims, 
              ylim = lims,
              xlab = "+SMARCA4 Log2FC", 
              ylab = "+SMARCA2 Log2FC")
dev.off()

pdf(file.path(".", "output", "figure3", "fig3b_wt_rescue_scatter_t910m.pdf"), width = 5, height = 5)
smoothScatter(log(atac_df$SMARCA4_FC+1, 2),
              log(atac_df$T910_FC+1, 2),
              colramp = colorRampPalette(RColorBrewer::brewer.pal(9, "Greys")),
              xlim = lims, 
              ylim = lims,
              xlab = "+SMARCA4 Log2FC", 
              ylab = "+SMARCA4 T910M Log2FC")
dev.off()

pdf(file.path(".", "output", "figure3", "fig3b_wt_rescue_scatter_k785r.pdf"), width = 5, height = 5)

smoothScatter(log(atac_df$SMARCA4_FC+1, 2), 
              log(atac_df$K785R_FC+1, 2),
              colramp = colorRampPalette(RColorBrewer::brewer.pal(9, "Greys")),
              xlim = lims, 
              ylim = lims,
              xlab = "+SMARCA4 Log2FC", 
              ylab = "+SMARCA4 K785R Log2FC")
dev.off()

