#09-RNA_seq_ChIP_seq_Integration

library(ProjectTemplate); load.project()
library(pheatmap)
source(file.path(".", "munge", "05-Gene_Expression.R"))
source(file.path(".", "munge", "06-Differential_Expression.R"))
source(file.path(".", "munge", "07-MACS_peaks.R"))
source(file.path(".", "munge", "08-Peak_to_Gene.R"))

### RNA-seq, ChIP-seq integration ###
color_map <- c("BAF activity-dependent" = "green",
               "BAF activity-independent" = "red",
               "PBAF activity-dependent" = "blue",
               "PBAF activity-independent" = "orange",
               "Residual" = "lightgray",
               "Neutral" = "white")

gexp_mat <- bin_deseq_complete %>%
  filter(Gene %in% c(bin_wt_up, bin_wt_down)) %>%
  select(Gene, Condition, log2FoldChange) %>%
  df.to.mat()

#Max value 12 for color bar
gexp_mat[gexp_mat > 12] <- 12

baf_annot <- rownames(gexp_mat) %>% 
  enframe(value = "Gene") %>% 
  left_join(baf_peak_to_gene) %>% 
  select(Peakset) %>% 
  replace_na(list(Peakset = "Neutral")) %>% 
  as.data.frame() %>% 
  set_rownames(rownames(gexp_mat))

pbaf_annot <- rownames(gexp_mat) %>% 
  enframe(value = "Gene") %>% 
  left_join(pbaf_peak_to_gene) %>% 
  select(Peakset) %>% 
  replace_na(list(Peakset = "Neutral")) %>% 
  as.data.frame() %>% 
  set_rownames(rownames(gexp_mat))

annot <- cbind(baf_annot %>% rename(BAF = Peakset), 
               pbaf_annot  %>% rename(PBAF = Peakset)) 
#Heatmaps
trees <- pheatmap(gexp_mat, 
                  show_rownames = F, 
                  annotation_row = annot, 
                  cluster_rows = T,
                  annotation_colors = list("BAF" = color_map,
                                           "PBAF" = color_map))

#output
pdf(file.path(".", "output", "figure6", "fig6a_gexp_heatmap.pdf"), height = 7, width = 5)
trees
dev.off()

#clusters
clusters <- cutree(trees$tree_row, k = 3) %>% 
  enframe(name = "Gene", value = "Cluster")

#Cluster identity boxplots
clusters  %>% 
  left_join(bin_deseq) %>% 
  separate(Condition, into = c("Cell_Line", "Condition")) %>% 
  mutate(Condition = factor(Condition, levels = c("K785R",
                                                  "T910M",
                                                  "SMARCA4",
                                                  "SMARCA2") )) %>% 
  filter(!is.na(Condition)) %>% 
  ggplot(aes(Condition, log2FoldChange)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~Cluster, nrow = 4) + 
  ylim(-10, 15) +
  geom_hline(yintercept = 1, linetype="dashed", color = "gray") +
  geom_hline(yintercept = -1, linetype="dashed", color = "gray") +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme_minimal() +
  ggsave(file.path(".", "output", "figure6", "fig6c_cluster_boxplot.pdf"), height = 5, width = 2.25)

#Genesets
baf_dep_genes <- baf_peak_to_gene %>%
  filter(Peakset == "BAF activity-dependent") %>% 
  pull(Gene)

baf_ind_genes <- baf_peak_to_gene %>%
  filter(Peakset == "BAF activity-independent") %>% 
  pull(Gene)

pbaf_dep_genes <- pbaf_peak_to_gene %>%
  filter(Peakset == "PBAF activity-dependent") %>% 
  pull(Gene)

pbaf_ind_genes <- pbaf_peak_to_gene %>%
  filter(Peakset == "PBAF activity-independent") %>% 
  pull(Gene)

target_genes <- list(baf_dep_genes,
                     baf_ind_genes,
                     pbaf_dep_genes,
                     pbaf_ind_genes)


#Enrichment tests
baf_cluster_counts <- clusters %>% 
  left_join(baf_peak_to_gene) %>% 
  group_by(Cluster) %>% 
  dplyr::count(Peakset) %>% 
  replace_na(list(Peakset = "Neutral")) %>% 
  mutate(Cluster_Total = sum(n)) %>% 
  spread(Peakset, n) %>% 
  gather(Peakset, Count, `BAF activity-dependent`, `BAF activity-independent`) %>% 
  ungroup()

pbaf_cluster_counts <- clusters %>% 
  left_join(pbaf_peak_to_gene) %>% 
  group_by(Cluster) %>% 
  dplyr::count(Peakset) %>% 
  replace_na(list(Peakset = "Neutral")) %>% 
  mutate(Cluster_Total = sum(n)) %>% 
  spread(Peakset, n) %>% 
  gather(Peakset, Count, `PBAF activity-dependent`, `PBAF activity-independent`) %>% 
  ungroup()

background_df <- tibble(Peakset = c("BAF activity-dependent",
                                    "BAF activity-independent",
                                    "PBAF activity-dependent",
                                    "PBAF activity-independent"),
                        Total = lengths(target_genes),
                        Background = length(bin_expressed) - Total)

fisher_df <- rbind(baf_cluster_counts, pbaf_cluster_counts) %>% 
  left_join(background_df) %>% 
  mutate(Fisher_Test = purrr::pmap(list(Count,
                                        Cluster_Total - Count,
                                        Total,
                                        Background), custom_fishers_greater)) %>% 
  mutate(Glances = map(Fisher_Test, broom::glance)) %>% 
  unnest(Glances) %>%
  group_by(Peakset) %>% 
  mutate(padj = p.adjust(p.value)) %>% 
  ungroup()

fisher_df  %>% 
  select(Cluster, Peakset, p.value, padj) %>% 
  mutate(Cluster = factor(Cluster, levels = c(1, 3, 2)) %>% 
           revalue(c("1" = "Downregulated",
                     "2" = "Upregulated",
                     "3" = "Highly upregulated"))) %>% 
  ggplot(aes(Peakset, -log(padj, 10), fill = Peakset)) +
  geom_col() +
  facet_wrap(~Cluster, ncol = 1) +
  coord_flip() +
  scale_fill_manual(values = c("BAF activity-dependent" = "#6ABD45",
                               "BAF activity-independent" = "#ED2024",
                               "PBAF activity-dependent" = "#3953A4",
                               "PBAF activity-independent" = "#FAA41A")) +
  BuenColors::pretty_plot() +
  theme(legend.position = "none") +
  ggsave(file.path(".", "output", "figure6", "fig6b_cluster_peak_enrichment.pdf"), height = 5, width = 3.5)

### Gorilla ###
cluster_genesets <- split(clusters$Gene, clusters$Cluster)

#Export genesets for GOrilla
cat(cluster_genesets[[1]], file= file.path(".", "output", "files", "1.tsv"), sep = "\n")
cat(cluster_genesets[[2]], file= file.path(".", "output", "files", "2.tsv"), sep = "\n")
cat(cluster_genesets[[3]], file= file.path(".", "output", "files", "3.tsv"), sep = "\n")

#Load cached results from GOrilla web portal
list(Downregulated = read_tsv(file.path(".", "data", "interim", "GOrilla", "1.xls")),
     Upregulated = read_tsv(file.path(".", "data", "interim", "GOrilla", "2.xls")),
     Highly_Upregulated = read_tsv(file.path(".", "data", "interim", "GOrilla", "3.xls"))) %>% 
  llply(function(x) x %>% slice(1:4)) %>% 
  enframe() %>% 
  unnest() %>% 
  ggplot(aes(rep(4:1, 3), -log(`FDR q-value`, 10), label = Description)) +
  geom_col() + 
  coord_flip() + 
  facet_wrap(~name, ncol =  1) + 
  geom_text() + 
  theme_minimal() +
  ggsave(file.path(".", "output", "figure6", "fig6c_cluster_go_enrich.pdf"), width = 2.25, height= 5)

