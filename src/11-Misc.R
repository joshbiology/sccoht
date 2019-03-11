#11-Misc.R

library(ProjectTemplate)
load.project()

#IF quant
process_IF <- function(files) {
  files %>% 
    adply(1, function(x) read_tsv(file.path("./data/raw/IF", x))) %>% 
    mutate(ID_1 = X1 %>% as.character,
           ID_2 = Area %>% as.character) %>% 
    unite(ID, ID_1, ID_2) %>% 
    mutate(Stain = word(Label, 1, sep = "\\.") %>% str_sub(1, 3),
           CTCF = IntDen - Area * Mean) %>% 
    filter(Area > 0.15) %>% 
    select(ID, Stain, CTCF) %>% 
    spread(Stain, CTCF) %>% 
    mutate(brg = brg/dap,
           ss1 = ss1/dap)
}

if_data <- list(list.files("./data/raw/IF/", pattern = "1_.*xls"),
                list.files("./data/raw/IF/", pattern = "2_.*xls")) %>% 
  ldply(process_IF) %>% 
  filter(brg > -4)

if_cor <- cor(if_data$brg, if_data$ss1, method = "spearman")

if_data %>% 
  ggplot(aes(brg, ss1)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(data = tibble(brg = 0.25, ss1 = 2.4),
            aes(size = 20, label = paste("cor = ", round(if_cor, digits = 4)))) +
  scale_x_continuous("SMARCA4 Intensity", c(0, 1, 2, 2.5), c(0, 2.5)) +
  scale_y_continuous("SS18 Intensity", c(0, 1, 2, 2.5), c(0, 2.5)) + 
  theme(legend.position="none") +
  ggsave(file.path(".", "output", "suppfigure2", "sfig2c_IF_scatter.pdf"), width = 3, height = 3)


#### Bar plot showing the increase in peaks
peakset <- c("BIN67_Control2_ARID2",
             "BIN67_SMARCA4_ARID2",
             "BIN67_Control2_SS18",
             "BIN67_SMARCA4_SS18",
             "BIN67_Control2_Brg",
             "BIN67_SMARCA4_Brg",
             "BIN67_Control2_155",
             "BIN67_SMARCA4_155")

load(file.path(".", "cache", "bin_chip_peaks.RData"))

bin_chip_peaks[peakset] %>% 
  ldply(length) %>% 
  set_colnames(c("Sample", "Peaks")) %>% 
  separate(Sample, c("Cell_Line", "Condition", "Epitope"), sep = "_") %>% 
  ggplot(aes(Epitope, Peaks, fill = Condition)) +
  geom_col(position = "dodge") +
  ggsave(file.path(".", "output", "suppfigure2", "sfig2h_barplot_peak_number.pdf"), width = 5, height = 3)