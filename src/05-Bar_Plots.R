#05-Bar_Plots.R

library(ProjectTemplate)
load.project()
source(file.path(".", "munge", "05-Gene_Expression.R"))
source(file.path(".", "munge", "07-MACS_peaks.R"))

#Figure 3C
tibble(Number = lengths(complex_k785r[c(2, 3, 7, 8)]),
       Group = c("Activity-independent", "Activity-dependent", "Activity-independent", "Activity-dependent"),
       Epitope = c("BAF", "BAF", "PBAF", "PBAF") %>% factor(level = c("PBAF", "BAF"))) %>% 
  group_by(Epitope) %>% 
  mutate(Total = sum(Number)) %>% 
  mutate(Frac = Number / Total) %>% 
  ggplot(aes(Epitope, Frac, fill = Group)) +
  geom_col() +
  theme(legend.position = "bottom") +
  ggsave(file.path(".", "output", "figure3", "fig3c_Frac_Activity_Dependent.pdf"), width = 2.5, height = 4)

#fishers
fisher.test(matrix(lengths(complex_k785r[c(2, 3, 7, 8)]), nrow = 2))

#Figure 4A and 5A
#Use the peak-gene assignment logic for PBAF
tmp <- map_dfr(list(Ctrl = complex_k785r[4], 
                    Mut = complex_k785r[c(4, 7)], 
                    WT = complex_k785r[c(4, 7, 8)]), 
               function(x) assign_peakset_to_geneset(x, tss_4kb),
               .id = "Condition") %>% 
  group_by(Condition, Peakset) %>% 
  dplyr::count() %>% 
  mutate(Frac = n/length(tss_4kb)) %>% 
  select(-n) %>% 
  ungroup()

tmp2 <- tmp %>% 
  group_by(Condition) %>% 
  summarize(Frac = 1-sum(Frac)) %>% 
  mutate(Peakset = "Unoccupied") %>% 
  select(1, 3, 2) %>% 
  ungroup()

bind_rows(tmp, tmp2) %>% 
  ggplot(aes(Condition, y = Frac, fill = Peakset)) + 
  geom_col(position = "stack") +
  scale_fill_manual(values = c("Ctrl_Residual_Reference_Peaks" = "#C49A6C",
                               "PBAF_WT_Mut_Shared" = "#18BDC2",
                               "PBAF_WT_Only" = "#F3766E",
                               "Unoccupied" = "gray"))  +
  theme(legend.position = "NA") + 
  ggsave(file.path(".", "output", "figure5", "fig5a_pbaf_tss_barplots.pdf"), width = 2.75, height = 4)

#Use the peak-enhancer assignment logic for PBAF
tmp3 <- map_dfr(list(Ctrl = complex_k785r[4], 
                     Mut = complex_k785r[c(4, 2)], 
                     WT = complex_k785r[c(4, 2, 3)]), 
                function(x) assign_peakset_to_ref_set(x, enh),
                .id = "Condition") %>% 
  group_by(Condition, Peakset) %>% 
  dplyr::count() %>% 
  mutate(Frac = n/length(enh)) %>% 
  select(-n) %>% 
  ungroup()

tmp4 <- tmp3 %>% 
  group_by(Condition) %>% 
  summarize(Frac = 1-sum(Frac)) %>% 
  mutate(Peakset = "Unoccupied") %>% 
  select(1, 3, 2) %>% 
  ungroup()

bind_rows(tmp3, tmp4) %>% 
  mutate(Peakset = factor(Peakset, levels = c("Ctrl_Residual_Reference_Peaks",
                                              "BAF_WT_Mut_Shared",
                                              "BAF_WT_Only",
                                              "Unoccupied"))) %>% 
  ggplot(aes(Condition, y = Frac, fill = Peakset)) + 
  geom_col(position = "stack") +
  scale_fill_manual(values = c("Ctrl_Residual_Reference_Peaks" = "#C49A6C",
                               "BAF_WT_Mut_Shared" = "#18BDC2",
                               "BAF_WT_Only" = "#F3766E",
                               "Unoccupied" = "gray")) +
  theme(legend.position = "NA")  + 
  ggsave(file.path(".", "output", "figure4", "fig4a_baf_enh_barplots.pdf"), width = 2.75, height = 4)
