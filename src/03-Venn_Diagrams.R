#03-Venn Diagrams

library(ProjectTemplate)
load.project()
source(file.path(".", "munge", "05-Gene_Expression.R"))
source(file.path(".", "munge", "07-MACS_peaks.R"))

#Calculate overlaps
ctrl_overlap <- findOverlapsOfPeaks(single_epitopes[["Ctrl_ARID2_SMARCC1_Overlap"]], single_epitopes[["Ctrl_DPF2_SMARCC1_Overlap"]], maxgap = 100)

pdf(file.path(".", "output", "figure2", "fig2c_Ctrl_DPF2_ARID2_Overlap.pdf"), width = 5, height = 5)
makeVennDiagram(ctrl_overlap)
dev.off()

wt_overlap <- findOverlapsOfPeaks(single_epitopes[["WT_ARID2_SMARCA4_Overlap"]], single_epitopes[["WT_SS18_SMARCA4_Overlap"]], maxgap = 100)

pdf(file.path(".", "output", "figure2", "fig2e_WT_SS18_ARID2_Overlap.pdf"), width = 5, height = 5)
makeVennDiagram(wt_overlap)
dev.off()

