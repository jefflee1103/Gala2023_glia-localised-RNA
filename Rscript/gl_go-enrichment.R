# Run GO enrichment analysis on the glial-localised transcripts 
# Filtered for ID interest

# ----- Environment

library(topGO)
library(tidyverse)
# need to load appropriate species database
library(org.Dm.eg.db)
library(qs)
library(rrvgo)
library(ggrepel)
library(furrr)
library(simplifyEnrichment)
library(patchwork)
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
ont_string <- c("BP", "MF", "CC") %>% set_names() 

se_goid <- ont_string %>%
  map(~{ 
    filter(gl_go, bonferroni < 0.01 & foldEnrichment > 2) %>%
      filter(ontology == .x) %>%
      pull(GO.ID)
  })
se_simmat <- future_map(se_goid, ~ GO_similarity(.x, db = "org.Dm.eg.db"))
future_map(se_simmat, select_cutoff) %>% wrap_plots()
pdf("~/Desktop/se-test_fc2.0_cutoff0.75.pdf", width = 8, height = 8)
se_simplify <- map(se_simmat, ~ simplifyGO(.x, control = list(cutoff = 0.75)))
dev.off()



# ----- Plotting

gl_go %>%
  mutate(put_label = if_else(
    log2(foldEnrichment) > 1.5 & bonferroni < 0.01,
    Term,
    NA_character_)) %>% 
  ggplot(aes(x = log2(foldEnrichment), y = -log(bonferroni), colour = ontology)) +
  geom_point(alpha = 0.5, stroke = 0, size =3) +
  geom_text_repel(aes(label = put_label), max.overlaps = 10,
    colour = "gray50", cex = 2) +
  labs(title = "GO enrichment: 1700 vs full genome",
    subtitle = "Half-volcano plot",
    x = "log2FoldChange", y = "Adjusted p-value (-log10)") +
  coord_cartesian(xlim = c(0, 3.5)) + 
  theme_light() +
  facet_wrap(~ ontology)

ggsave("./output/graphics/go-enrichment_volcano.pdf", width = 20, height = 14)



