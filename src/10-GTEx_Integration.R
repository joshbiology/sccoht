#10-GTEx_Integration

library(ProjectTemplate)
load.project()
source(file.path(".", "munge", "05-Gene_Expression.R"))
source(file.path(".", "munge", "06-Differential_Expression.R"))


### Integrating GTEX gene expression with BIN data ###

#Figure 6D
read_tsv(file.path(".", "data", "interim", "gtex", "GTex_v_Tumor_rank_table.txt")) %>% 
  mutate(Rank = nrow(.):1) %>% 
  ggplot(aes(Rank, -diff)) +
  geom_line() +
  ggsave(file.path(".", "output", "figure6", "fig6d_diff_plot.pdf"), width = 4, height = 6)

#Figure 6E
ovary_gtex <- map(list.files(file.path(".", "data", "interim", "gtex"),
                             pattern = "genes_in_tumor.txt",
                             full.names = T),
                  function(x)read_tsv(x, col_names = F)) %>% 
  set_names(c("Up", "Down")) %>% 
  enframe() %>% 
  unnest() %>% 
  set_colnames(c("Direction", "Gene"))

ovary_gtex %>% dplyr::count(Direction)

ovary_gtex %>%
  mutate(Direction = factor(Direction, levels = c("Up", "Down"))) %>% 
  left_join(bin_deseq) %>% 
  mutate(Condition = case_when(Condition == "SCCOHT_SMARCA4" ~"SMARCA4 WT Rescue",
                               Condition == "SCCOHT_T910M" ~"SMARCA4 T910M Rescue",
                               Condition == "SCCOHT_SMARCA2" ~"SMARCA2 WT Rescue",
                               TRUE ~ "SMARCA4 K785R Rescue") %>% 
           factor(levels = c("SMARCA4 K785R Rescue", "SMARCA4 T910M Rescue", "SMARCA4 WT Rescue", "SMARCA2 WT Rescue"))) %>% 
  filter(Condition %in% c("SMARCA4 WT Rescue", "SMARCA2 WT Rescue")) %>% 
  ggplot(aes(log2FoldChange, -log(pvalue, 10), color = Direction)) +
  geom_point(size= 0.5) +
  xlim(-5, 5) +
  ylim(0, 100) +
  scale_color_manual(values = c(Up = "#ED1C24", Down = "#3853A4")) +
  facet_grid(Direction~Condition) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Log2 Fold Change RPKM",
       y = "Significance (-Log10 P-val)") + 
  geom_vline(xintercept = c(-1, 1), color = "gray", linetype = 2) +
  ggsave(file.path(".", "output", "figure6", "fig6e_tumor_gexp_updated.pdf"), width = 5, height = 5)

#Figure 6F
ovary_gtex %>%
  mutate(Direction = factor(Direction, levels = c("Up", "Down"))) %>% 
  inner_join(bin_deseq %>% filter(Condition == "SCCOHT_SMARCA4")) %>% 
  filter(abs(log2FoldChange) > 1) %>% 
  mutate(Direction_BIN = log2FoldChange > 0)%>% 
  dplyr::count(Direction, Direction_BIN)

fisher.test(matrix(c(311, 56, 58, 292), nrow = 2), alternative = "g")
