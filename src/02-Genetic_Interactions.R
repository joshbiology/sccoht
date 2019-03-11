#02_Genetic_Interactions.R
#Original author: Benedikt Rauscher (MINGLE, Rauscher et al. 2018 MSB)
#Modified to focus on Dual ATPase loss

library(ProjectTemplate)
load.project()
source(file.path(".", "munge", "02-Dependency.R"))
source(file.path(".", "munge", "03-Genesets.R"))

library(ggsignif)
library(ggrepel)

#Subunit loss table: Robin Meyers, from Pan, Meyers et al. Cell Systems 2018
subunit_loss_table <- read_tsv("./data/interim/loss_calls/SWISNF_Subunit_deficiency.tsv")

# Create BAF subunit loss table as matrix
somatic_filtered <- subunit_loss_table %>% 
  select(Query = Subunit, Cell_Line = Sample, Status = Class1) %>% 
  spread(Query, Status) %>% 
  as.data.frame() %>% 
  set_rownames(.$Cell_Line) %>% 
  select(-Cell_Line) %>% 
  as.matrix()

dual_annot <- cbind(Dual_ATPase = somatic_filtered[,"SMARCA4"] * somatic_filtered[,"SMARCA2"],
                    Dual_ARID1 = somatic_filtered[,"ARID1A"] * somatic_filtered[,"ARID1B"])

#Create dual loss annotations w/ literature
somatic_filtered[,"SMARCA4"][dual_annot[,"Dual_ATPase"] == 1] <- 0
somatic_filtered[,"SMARCA2"][dual_annot[,"Dual_ATPase"] == 1] <- 0

somatic_filtered[,"ARID1A"][dual_annot[,"Dual_ARID1"] == 1] <- 0
somatic_filtered[,"ARID1B"][dual_annot[,"Dual_ARID1"] == 1] <- 0

somatic_filtered["MIAPACA2_PANCREAS","ARID1B"] <- 0

somatic_filtered <- somatic_filtered %>% 
  cbind(dual_annot) %>%
  as.data.frame() %>% 
  mutate(Cell_Line = rownames(.)) %>% 
  gather(Query, Status, 1:(ncol(.)-1)) %>% 
  filter(Status == 1)

#Inspect for sensical outputs
somatic_filtered %>% 
  filter(Cell_Line == "MIAPACA2_PANCREAS")

#Manually add additional literature curated loss events
manual_annot <- tibble(Cell_Line = c("MIAPACA2_PANCREAS", "CHLA266_SOFT_TISSUE", "COGAR359_SOFT_TISSUE", "BIN67_OVARY"),
                       Query = c("Dual_ARID1", "SMARCB1", "SMARCB1", "Dual_ATPase"),
                       Status = c(1, 1, 1, 1))

somatic_filtered <- somatic_filtered %>% 
  rbind(manual_annot)

colnames(somatic_filtered) <- c("cellline", "symbol", "status")

somatic_filtered <- somatic_filtered %>% 
  filter(symbol %in% c(swisnf_genes, "Dual_ARID1", "Dual_ATPase"),
         cellline %in% colnames(avana_dep))

#FOCUSED ANALYSES ON Dual loss of ATPase
mutated_BAF_cell_lines <- somatic_filtered %>% 
  filter(symbol != "Dual_ATPase") %>% 
  pull(cellline) %>% 
  unique()

avana_df <- avana_dep[avana_analysis_genes,] %>%
  t() %>%
  as_tibble() %>%
  tibble::add_column(Cell_Line = colnames(avana_dep), .before = 1) %>%
  tibble::add_column(Lineage = word(colnames(avana_dep), sep = "_", 2, -1), .after = 1) %>% 
  gather(Gene, Dep, 3:ncol(.))

#Curate list to include only dual_ATPase
somatic_filtered <- somatic_filtered %>% 
  filter(symbol == "Dual_ATPase") %>% 
  left_join(avana_df %>% filter(Gene == "SMARCA4") %>% select(cellline = Cell_Line, Dep)) %>% 
  filter(Dep > -0.5) %>% 
  select(-Dep)


#Boxplots / Jitter plots of differential essentiality
gns <- c("SMARCB1", "SMARCE1", "ARID1A", "PBRM1", "ARID2", "SMARCA4", "EED", "SUZ12", "EZH2")
swisnf_dep <- avana_dep  %>% 
  mat.to.df(row.name = "Gene", col.name = "Cell_Line", dat.name = "Dep") %>%
  filter(Gene %in% gns) %>% 
  left_join(somatic_filtered, by = c("Cell_Line" = "cellline")) %>% 
  mutate(Gene = factor(Gene, levels = gns))

swisnf_dep %>% 
  ggplot(aes(Gene, Dep, fill = symbol)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggsave(file.path(".", "output", "suppfigure1", "fig1h_ess_boxplots.pdf"), width = 6, height = 4)

#Create bfs_norm
bfs_norm <- avana_dep %>%
  mat.to.df("Gene", "Cell_Line", "Dep") %>% 
  filter(Gene %in% avana_analysis_genes,
         !(Cell_Line %in% mutated_BAF_cell_lines)) %>% 
  select(symbol = Gene, cellline = Cell_Line, BF = Dep)

#### Create combinations
combis_batch <- expand.grid(bfs_norm %>% .$symbol %>% unique,
                            unique(somatic_filtered$symbol)) %>% as_tibble() %>%
  dplyr::rename(target=Var1, query=Var2) %>%
  mutate(target=as.character(target), query=as.character(query)) %>% 
  filter(query == "Dual_ATPase")

## merge all data into big data frame for testing.
combis_bf_batch <- combis_batch %>%
  inner_join(bfs_norm %>% dplyr::rename(target=symbol)) %>%
  left_join(somatic_filtered %>% mutate(mutated=1) %>%
              dplyr::select(symbol, cellline, mutated) %>% dplyr::rename(query=symbol)) %>%
  mutate(type=ifelse(is.na(mutated), 'target', 'combi')) %>%
  dplyr::select(-c(mutated)) %>% distinct()

### Calculate statistic

sgi_filtered <- combis_bf_batch %>% 
  group_by(target, query) %>% 
  do(broom::tidy(wilcox.test(.$BF ~ .$type))) %>% 
  ungroup() %>% 
  mutate(q.value = p.adjust(.$p.value, method = "fdr"))

##### Calculate pi-scores
#Process dataframe into: Target, Query (Gene lost in particular cell line), Dependency
### Dependency dataframes. Tidy, only analysis genes included


avana_query <- avana_df %>%
  dplyr::rename(Target = Gene) %>% 
  #filter(Target %in% union(swisnf_genes, prc2_genes)) %>%
  left_join(somatic_filtered, by = c("Cell_Line" = "cellline")) %>% 
  dplyr::rename(cellline = Cell_Line, target = Target, query = symbol, BF = Dep)

#Create pi scores
pl_out <- avana_query %>% 
  dplyr::group_by(target, query) %>% 
  dplyr::summarise(fitness=mean(BF)) %>% 
  ungroup() %>% 
  spread(query,fitness) %>%
  data.frame() %>% 
  `rownames<-`(.$target) %>% 
  dplyr::select(-target) %>%
  as.matrix() %>% 
  HD2013SGI::HD2013SGImaineffects()

all_effects <- pl_out$pi %>% melt %>%
  `colnames<-`(c('target', 'query', 'pi')) %>% tbl_df %>%
  inner_join(pl_out$targetMainEffect %>% melt %>% mutate(target=rownames(.)) %>%
               dplyr::rename(target_main=value)) %>%
  inner_join(cbind(pl_out$queryMainEffect, pl_out$pi %>% colnames()) %>% tbl_df %>%
               `colnames<-`(c('query_main', 'query')) %>%
               mutate(query_main=as.numeric(query_main)))


pi_scores <- avana_query %>% dplyr::select(target, query, BF) %>%
  mutate(fitness=BF) %>% dplyr::select(-BF) %>%
  group_by(target, query) %>% 
  dplyr::summarise(fitness=mean(fitness)) %>% 
  ungroup() %>% 
  mutate(fitness=scale(fitness)[,1]) %>%
  dplyr::rename(measured=fitness) %>%
  inner_join(all_effects)

#### Plotting and annnotations
pi_threshold <- 0.15

all_interactions <- pi_scores %>% 
  dplyr::inner_join(sgi_filtered) %>%
  dplyr::inner_join(
    avana_query %>% group_by(target, query) %>%
      dplyr::summarise(var_bf=var(BF)) %>% ungroup() %>% distinct()
  )


#all_interactions %>% 
#	write_tsv("./output/files/CRISPR_dual_ATPase_epistasis.tsv")

genes_of_interest_df <- tibble(target = c("SMARCB1", "SMARCE1", "ARID1A", "PBRM1", "ARID2", "SMARCA4", "SMARCA2", "EED", "SUZ12", "EZH2", "SMARCD1"),
                               Geneset = c("BAF", "BAF", "BAF", "PBAF", "PBAF", "Lost", "Lost", "PRC2", "PRC2", "PRC2", "Shared"))
dual_atpase_inter <- all_interactions %>% dplyr::filter(query == "Dual_ATPase")

plt_df <- dual_atpase_inter %>% 
  mutate(
    type=ifelse(q.value < 0.25 & pi > 0, 'positive',
                ifelse(q.value < 0.25 & pi < 0, 'negative', 'normal')),
    type=factor(type, levels=c('negative', 'normal', 'positive'))) %>% 
  left_join(genes_of_interest_df)

#ignore for now
plt_df %>% 
  ggplot(aes(x=-log10(p.value), y=pi)) + 
  geom_point(aes(colour=type)) +
  scale_color_manual(values = c(negative = "blue", normal = "gray", "positive" = "yellow")) +
  geom_point(data = subset(plt_df, !is.na(Geneset)), aes(x=-log10(p.value), y=pi, size=2)) +
  geom_hline(yintercept=0, colour='#444444') + theme_classic() +
  geom_text_repel(data = subset(plt_df, !is.na(Geneset)), size=3, aes(label=target)) +
  theme(legend.position='bottom') + 
  xlab('p-value [-log10]') + 
  ylab('Pi-score') + 
  ggtitle("Dual_ATPase") +
  ggsave(file.path(".", "output", "figure1", "fig1f_dual_loss_ATPase_CRISPR.pdf"), width = 6, height = 5)

