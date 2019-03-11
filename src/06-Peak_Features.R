#06-Peak_Features.R
library(ProjectTemplate)
load.project()
source(file.path(".", "munge", "05-Gene_Expression.R"))
source(file.path(".", "munge", "07-MACS_peaks.R"))

#Figure 3F

library(matrixStats)


### Read-count boxplots ###

# Figure 4G: Boxplots for histones

read_tsv(file.path(".", "data", "interim", "rsubread", "hist_counts.tsv")) %>% 
  filter(Epitope != "ATAC") %>% 
  ggplot(aes(x = Epitope, y = Dat, fill = Dataset)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(0,4) +
  facet_wrap(~Condition, ncol = 2) +
  labs(y = "RPKM", x = "Enhancer Mark") +
  scale_fill_manual(values = c("Activity-independent" = "#00BFC4", 
                               "Activity-dependent" = "#F8766D")) +
  coord_flip() +
  theme_minimal() +
  ggsave(file.path(".", "output", "figure4", "fig4c_histone_marks_boxplot.pdf"), width = 6, height =4)

#Figure 4C

read_tsv(file.path(".", "data", "interim", "rsubread", "ss18_counts.tsv")) %>%
  mutate(.id = factor(.id, levels = c("Activity-independent", "Activity-dependent")),
         Assay = factor(Assay, levels = c("ChIP", "ATAC")),
         Condition = factor(Condition, levels = c("Control", "SMARCA4", "T910M", "K785R"))) %>% 
  ggplot(aes(x = Condition, y = Dat, fill = `.id`)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(0,4.5) +
  facet_wrap(.id ~ Assay, nrow = 2) + 
  theme_minimal() +
  labs(y = "RPKM", x = "Condition") +
  scale_fill_manual(values = c("#18BDC2", "#F3766E")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  ggsave(file.path(".", "output", "figure3", "fig3f_activity_readcounts.pdf"), width = 4, height = 5)

#Figure 5C
library(ggridges)

#Make density plot for ARID2
bin_broad_peaks[c("BIN67_K785R_ARID2_ZM101",
                  "BIN67_Control_ARID2_ZM087",
                  "BIN67_SMARCA4_ARID2_BIN9",
                  "BIN67_T910M_ARID2_ZM090")] %>% 
  magrittr::set_names(c("K785R", "Ctrl", "WT", "T910M")) %>% 
  llply(width) %>% 
  ldply(data.frame) %>% 
  set_colnames(c("Dataset", "Width")) %>% 
  mutate(Dataset = factor(Dataset, levels = c("Ctrl",
                                              "K785R",
                                              "T910M",
                                              "WT"))) %>% 
  ggplot(aes(Width, Dataset, fill = Dataset)) +
  geom_density_ridges_gradient(scale = 1) +
  scale_x_log10(limits = c(100, 10000)) + 
  coord_flip() + 
  scale_fill_manual(values = c(K785R = "#F0F9E8",
                               T910M = "#CCEBC5",
                               WT = "#78b5e8",
                               Ctrl = "#A9A9A9")) + 
  theme(legend.position = "none") +
  ggsave(file.path(".", "output", "figure5", "fig5c_arid2_peak_width_ridges.pdf"), height = 4, width = 3)

