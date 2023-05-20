# Run GO enrichment analysis on the glial-localised transcripts 
# Filtered for ID interest

# ----- Environment

# library(topGO)
library(tidyverse)
# need to load appropriate species database
library(org.Dm.eg.db)
library(qs)
# library(rrvgo)
library(ggrepel)
library(furrr)
library(simplifyEnrichment)
library(patchwork)
library(ggview)
plan(multisession, workers = 4)

source("./Rscript/runTopGO.R")

# ----- Get gene vectors

## Target genes
gl_nd <- qread("./output/glia-localised-neurodegenerative-disease-incorporated.qs")
id_interest <- gl_nd %>%
  filter(rna_in_protrusion >= 8) %>%
  filter(`subperineurial glial cell` == TRUE | 
           `perineurial glial sheath` == TRUE |
           `adult brain perineurial glial cell` == TRUE |
           `ensheathing glial cell` == TRUE)

targetGenes <- id_interest$dmel_gene_id

## Background genes
ens99 <- read_csv("./data/Flybase/Dmel_tx2gene_ENSEMBL_v99.csv") %>%
  dplyr::select(gene_id, gene_name) %>% distinct()

backgroundGenes <- ens99$gene_id

# ----- Run analysis

## Topgo
gl_go <- runTopGO_comprehensive("fly", backgroundGenes, targetGenes)
View(gl_go)

write_csv(gl_go, "./output/analysis/gl_go-enrichment.csv")

# ----- Reduce GO terms 

## Calculate similarity matrix 
gl_go <- read_csv("./output/analysis/gl_go-enrichment.csv")
ont_string <- c("BP", "MF", "CC") %>% set_names() 

rrvgo_filt_foldenrichment <- ont_string %>%
  future_imap(~ {
    go_df <- filter(gl_go, ontology == .y) %>%
      filter(foldEnrichment > 2 & bonferroni < 0.01)
    go_id_vector <- go_df$GO.ID
    scores <- setNames(go_df$foldEnrichment, go_df$GO.ID)
    simMatrix <- calculateSimMatrix(go_id_vector, orgdb = "org.Dm.eg.db", ont = .y, method = "Rel")
    reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold = 0.75, orgdb = "org.Dm.eg.db")
    output <- list(simMatrix, reducedTerms) %>% set_names(c("simMatrix", "reducedTerms"))
    return(output)
  })

rrvgo_filt_bonferroni <- ont_string %>%
  future_imap(~ {
    go_df <- filter(gl_go, ontology == .y) %>%
      filter(foldEnrichment > 2 & bonferroni < 0.01)
    go_id_vector <- go_df$GO.ID
    scores <- setNames(-log10(go_df$bonferroni), go_df$GO.ID)
    simMatrix <- calculateSimMatrix(go_id_vector, orgdb = "org.Dm.eg.db", ont = .y, method = "Rel")
    reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold = 0.75, orgdb = "org.Dm.eg.db")
    output <- list(simMatrix, reducedTerms) %>% set_names(c("simMatrix", "reducedTerms"))
    return(output)
  })

scatterPlot(rrvgo_filt_foldenrichment$BP$simMatrix, rrvgo_filt_foldenrichment$BP$reducedTerms)
treemapPlot(rrvgo_filt_foldenrichment$BP$reducedTerms)
treemapPlot(rrvgo_filt_bonferroni$BP$reducedTerms)

## Try simplify enrichment
gl_go <- read_csv("./output/analysis/gl_go-enrichment.csv")
se_cutoffs <- list(c(0.01, 2), c(0.01, 1.5), c(0.01, 2)) %>%
  set_names(c("BP", "MF", "CC"))

ont_string <- c("BP", "MF", "CC") %>% set_names() 

se_goid <- ont_string %>%
  map(~{ 
    filter(gl_go, ontology == .x) %>%
    filter(bonferroni < (pluck(se_cutoffs, .x))[1] & 
             foldEnrichment > (pluck(se_cutoffs, .x))[2]) %>%
      
      pull(GO.ID)
  })

se_simmat <- future_map(se_goid, ~ GO_similarity(.x, db = "org.Dm.eg.db"))
future_map(se_simmat, select_cutoff) %>% wrap_plots()

## BP
bp_clusters <- simplifyGO(se_simmat$BP, control = list(cutoff = 0.75), plot = FALSE)
plot <- ht_clusters(se_simmat$BP, bp_clusters$cluster, 
            min_term = 10,
            max_words = 6,
            order_by_size = TRUE, 
            bg_gp = gpar(fill = "gray92", col = "#AAAAAA"),
            word_cloud_grob_param = list(max_width = unit(40, "mm")), 
            fontsize_range = c(2, 9),
            col = sequential_hcl(n = 30, palette = "Reds 3", rev = TRUE),
            exclude_words = c("cellular", "involved", "process", "regulation", "cell", "metabolic", "transport", "positive", "negative"))
pdf("/Users/sjoh4548/OneDrive - Nexus365/JL_upload/bp.pdf", width = 13 * 0.3937, height = 6.3 * 0.3937)
plot
dev.off()

bp_clusters_summary <- bp_clusters %>%
  mutate(id = str_replace(id, ":", "")) %>%
  group_by(cluster) %>%
  summarise(go_count = n(),
    go_ids = list(id)) %>%
  filter(go_count >= 10) %>%
  arrange(desc(go_count)) 

bp_clusters_go_ids <- bp_clusters_summary$go_ids %>%
  set_names(bp_clusters_summary$cluster)

bp_clusters_go_ids %>%
  map(function(clusters){
    clusters %>%
      map(function(id){
        df <- id_interest %>%
          filter(str_detect(go_biological_process, id))
        if(nrow(df) == 0){
          NA_character_
        } else {
          pull(df, dmel_gene_id)
        }
      }) %>% 
      purrr::reduce(append) %>%
      unique() %>%
      .[!is.na(.)]
  }) %>%
  map(~ length(.x))

## MF
mf_clusters <- simplifyGO(se_simmat$MF, control = list(cutoff = 0.9), plot = FALSE)
plot <- ht_clusters(se_simmat$MF, mf_clusters$cluster, 
                    min_term = 5,
                    max_words = 6,
                    order_by_size = TRUE, 
                    bg_gp = gpar(fill = "gray92", col = "#AAAAAA"),
                    word_cloud_grob_param = list(max_width = unit(40, "mm")), 
                    fontsize_range = c(4, 9),
                    col = sequential_hcl(n = 30, palette = "Greens 3", rev = TRUE),
                    exclude_words = c("adenyl", "guanyl", "activity", "small", "molecule"),
                    )
pdf("/Users/sjoh4548/OneDrive - Nexus365/JL_upload/mf.pdf", width = 13 * 0.3937, height = 6.3 * 0.3937)
plot
dev.off()

mf_clusters_summary <- mf_clusters %>%
  mutate(id = str_replace(id, ":", "")) %>%
  group_by(cluster) %>%
  summarise(go_count = n(),
    go_ids = list(id)) %>%
  filter(go_count >= 5) %>%
  arrange(desc(go_count)) 

mf_clusters_go_ids <- mf_clusters_summary$go_ids %>%
  set_names(mf_clusters_summary$cluster)

mf_clusters_go_ids %>%
  map(function(clusters){
    clusters %>%
      map(function(id){
        df <- id_interest %>%
          filter(str_detect(go_molecular_function, id))
        if(nrow(df) == 0){
          NA_character_
        } else {
          pull(df, dmel_gene_id)
        }
      }) %>% 
      purrr::reduce(append) %>%
      unique() %>%
      .[!is.na(.)]
  }) %>%
  map(~ length(.x))

## CC
cc_clusters <- simplifyGO(se_simmat$CC, control = list(cutoff = 0.9), plot = FALSE)
plot <- ht_clusters(se_simmat$CC, cc_clusters$cluster, 
                    min_term = 10,
                    max_words = 6,
                    order_by_size = TRUE, 
                    bg_gp = gpar(fill = "gray92", col = "#AAAAAA"),
                    word_cloud_grob_param = list(max_width = unit(40, "mm")), 
                    fontsize_range = c(4, 9),
                    col = sequential_hcl(n = 30, palette = "Blues 3", rev = TRUE),
                    exclude_words = c("mitochondrial", "large", "small", "complex"),
)
pdf("/Users/sjoh4548/OneDrive - Nexus365/JL_upload/cc.pdf", width = 13 * 0.3937, height = 6.3 * 0.3937)
plot
dev.off()

cc_clusters_summary <- cc_clusters %>%
  mutate(id = str_replace(id, ":", "")) %>%
  group_by(cluster) %>%
  summarise(go_count = n(),
    go_ids = list(id)) %>%
  filter(go_count >= 10) %>%
  arrange(desc(go_count)) 

cc_clusters_go_ids <- cc_clusters_summary$go_ids %>%
  set_names(cc_clusters_summary$cluster)

cc_clusters_go_ids %>%
  map(function(clusters){
    clusters %>%
      map(function(id){
        df <- id_interest %>%
          filter(str_detect(go_cellular_component, id))
        if(nrow(df) == 0){
          NA_character_
        } else {
          pull(df, dmel_gene_id)
        }
      }) %>% 
      purrr::reduce(append) %>%
      unique() %>%
      .[!is.na(.)]
  }) %>%
  map(~ length(.x))



##
gl <- read_tsv("/Users/sjoh4548/Desktop/glia-localised-rna.txt")
id_interest <- gl %>%
  filter(rna_in_protrusion >= 8) %>%
  filter(`ensheathing glial cell` == TRUE | `perineurial glial sheath` == TRUE | `subperineurial glial cell` == TRUE | `adult brain perineurial glial cell` == TRUE)

id_interest %>% colnames()

bp_clusters %>%
  group_by(cluster) %>%
  summarise(id)



# ----- Plotting

gl_go %>%
  mutate(ontology = fct_relevel(ontology, c("BP", "MF", "CC"))) %>%
  mutate(put_label = if_else(
    log2(foldEnrichment) > 1.5 & bonferroni < 0.01,
    Term,
    NA_character_)) %>% 
  ggplot(aes(x = log2(foldEnrichment), y = -log(bonferroni), colour = ontology)) +
  geom_point(alpha = 0.5, stroke = 0, size = 1.5) +
  geom_text_repel(aes(label = put_label), max.overlaps = 10,
    colour = "gray50", cex = 1.5) +
  labs(title = "Gene Ontology enrichment: 1700 glia-localised transcripts vs full genome",
    subtitle = "Half-volcano plot - Top 20 (foldchange) terms are labelled",
    x = "log2FoldChange", y = "Adjusted p-value (-log10)", 
    colour = "GO ontology") +
  coord_cartesian(xlim = c(0, 3.5)) + 
  theme_light(base_size = 6) +
  facet_wrap(~ ontology)

ggview(width = 18, height = 15, units = "cm")

ggsave("./output/graphics/go-enrichment_volcano.pdf", 
       width = 18, height = 12, units = "cm",
       useDingbats = FALSE)



######################################################################
# With clusterprofiler
######################################################################
library(tidyverse)
library(qs)
library(patchwork)
library(org.Dm.eg.db)
library(furrr)
library(rstatix)
library(colorspace)
plan(multisession, workers = 4)

library(clusterProfiler)
library(enrichplot)


# Get gene lists

## Target genes
gl_nd <- qread("./output/glia-localised-neurodegenerative-disease-incorporated.qs")
id_interest <- gl_nd %>%
  filter(rna_in_protrusion >= 8) %>%
  filter(`subperineurial glial cell` == TRUE | 
           `perineurial glial sheath` == TRUE |
           `adult brain perineurial glial cell` == TRUE |
           `ensheathing glial cell` == TRUE)

targetGenes <- id_interest$dmel_gene_id

## Background genes
ens99 <- read_csv("./data/Flybase/Dmel_tx2gene_ENSEMBL_v99.csv") %>%
  dplyr::select(gene_id, gene_name) %>% distinct()

backgroundGenes_ens99 <- ens99$gene_id

backgroundGenes_mmus_convertible <- gl_nd$dmel_gene_id %>% unique()

# Cluster profilier
ontology <- c("BP", "MF", "CC") %>% set_names()

enrichgo_ens99 <- ontology %>%
  map(~ enrichGO(
    gene = targetGenes,
    universe = backgroundGenes_ens99,
    OrgDb = org.Dm.eg.db,
    keyType = "ENSEMBL",
    ont = .x,
    pvalueCutoff = 0.01,
    readable = TRUE
  ) %>% pairwise_termsim()) 

enrichgo_mmusdmel <- ontology %>%
  map(~ enrichGO(
    gene = targetGenes,
    universe = backgroundGenes_mmus_convertible,
    OrgDb = org.Dm.eg.db,
    keyType = "ENSEMBL",
    ont = .x,
    pvalueCutoff = 0.01,
    readable = TRUE
  ) %>% pairwise_termsim()) 


enrichgo_ens99$BP@result %>% View()

enrichgo_ens99$BP@result %>%
  filter(Description == "autophagy") %>%
  pull(geneID) %>%
  str_split_1(pattern = "/") %>% sort()

enrichgo_ens99$BP %>%
  emapplot()

cnetplot(enrichgo_ens99$MF, showCategory = 15, cex_label_category = 1.5, cex_category = 2, cex_label_gene = 0.8, layout = "kk")
ggsave("~/Desktop/1700_MF_ens99.pdf", height = 40, width = 40)

cnetplot(enrichgo_mmusdmel$BP, showCategory = 20, cex_label_category = 1.5, cex_category = 2, cex_label_gene = 0.8, layout = "kk")
ggsave("~/Desktop/1700_BP_mmusdmel.pdf", height = 40, width = 40)

goplot(enrichgo_mmusdmel$BP)


goplot(ego)

