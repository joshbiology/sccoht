#08-RNA_seq_Analysis

library(ProjectTemplate)
load.project()
source(file.path(".", "munge", "05-Gene_Expression.R"))
source(file.path(".", "munge", "06-Differential_Expression.R"))

### Volcano plots ###
fdr_threshold <- 0.05
FC_threshold <- 2
bin_deseq %>% 
  filter(Gene %in% bin_expressed) %>% 
  mutate(Status = case_when(padj < fdr_threshold & log2FoldChange > FC_threshold ~ "Upregulated",
                            padj < fdr_threshold  & log2FoldChange < -FC_threshold ~ "Downregulated",
                            TRUE ~ "Not Significant")) %>% 
  mutate(Condition = case_when(Condition == "SCCOHT_SMARCA4" ~"SMARCA4 WT Rescue",
                               Condition == "SCCOHT_T910M" ~"SMARCA4 T910M Rescue",
                               Condition == "SCCOHT_SMARCA2" ~"SMARCA2 WT Rescue",
                               TRUE ~ "SMARCA4 K785R Rescue") %>% 
           factor(levels = c("SMARCA4 K785R Rescue", "SMARCA4 T910M Rescue", "SMARCA4 WT Rescue", "SMARCA2 WT Rescue"))) %>% 
  ggplot(aes(log2FoldChange, -log(pvalue, 10), color = Status)) +
  geom_point(size= 0.5) +
  xlim(-10, 10) +
  ylim(0, 150) +
  scale_color_manual(values = c(Upregulated = "#ED1C24", Downregulated = "#3853A4", `Not Significant` = "#B0B0B1")) +
  facet_wrap(~Condition, nrow = 1) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Log2 Fold Change RPKM",
       y = "Significance (-Log10 P-val)") +
  ggsave(file.path(".", "output", "suppfigure5", "sfig5a_gexp_volcano.png"), width = 7, height = 3, dpi = 300)

#Gene numbers
bin_deseq %>% 
  filter(Gene %in% bin_expressed) %>% 
  mutate(Status = case_when(padj < fdr_threshold & log2FoldChange > FC_threshold ~ "Upregulated",
                            padj < fdr_threshold  & log2FoldChange < -FC_threshold ~ "Downregulated",
                            TRUE ~ "Not Significant")) %>% 
  group_by(Condition, Status) %>% 
  dplyr::count()



