#08-Peak_to_Gene

gr1 <- complex_k785r[["BAF_WT_Only"]]
gr2 <- complex_k785r[["BAF_WT_Mut_Shared"]]
gr3 <- complex_k785r[["Ctrl_Residual_Reference_Peaks"]]
gr4 <- complex_k785r[["PBAF_WT_Only"]]
gr5 <- complex_k785r[["PBAF_WT_Mut_Shared"]]

#Peak to Gene mapping
#BAF
baf_peak_to_gene <- list("BAF activity-dependent" = gr1, 
                         "BAF activity-independent" = gr2) %>% 
  assign_peakset_to_geneset(tss_30kb)

pbaf_peak_to_gene <- list("PBAF activity-dependent" = gr4, 
                          "PBAF activity-independent" = gr5)  %>% 
  assign_peakset_to_geneset(tss_30kb)

baf_ind_targets <- baf_peak_to_gene %>% 
  filter(Peakset == "BAF activity-independent") %>% 
  pull(Gene)

baf_dep_targets <- baf_peak_to_gene %>% 
  filter(Peakset == "BAF activity-dependent") %>% 
  pull(Gene)

pbaf_ind_targets <- pbaf_peak_to_gene %>% 
  filter(Peakset == "PBAF activity-independent") %>% 
  pull(Gene)

pbaf_dep_targets <- pbaf_peak_to_gene %>% 
  filter(Peakset == "PBAF activity-dependent") %>% 
  pull(Gene)
