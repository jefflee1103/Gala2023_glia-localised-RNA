library(tidyverse)
library(qs)
library(biomaRt)

## Get 1700 genes
gl_nd <- qread("./output/glia-localised-neurodegenerative-disease-incorporated.qs")
id_interest <- gl_nd %>%
  filter(rna_in_protrusion >= 8) %>%
  filter(`subperineurial glial cell` == TRUE | 
           `perineurial glial sheath` == TRUE |
           `adult brain perineurial glial cell` == TRUE |
           `ensheathing glial cell` == TRUE)

## Get SFARI genes 
sfari_raw <- read_csv("./data/SFARI/SFARI-Gene_genes_11-07-2022release_11-15-2022export.csv") %>%
  janitor::clean_names()

## Fill in missing human ensembl ids
hsap_tx2gene <- read_csv("./data/SFARI/Hsap_tx2gene_ENSEMBL_v99.csv") %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct()

missing_ids <- sfari_raw %>%
  filter(is.na(ensembl_id)) %>%
  pull(gene_symbol)

missing_id_matches <- hsap_tx2gene %>%
  filter(gene_name %in% missing_ids) %>%
  filter(str_detect(gene_id, "ENSG")) %>%
  group_by(gene_name) %>%
  slice_head(n = 1)

sfari_cleaned <- bind_rows(
  filter(sfari_raw, !is.na(ensembl_id)),
  filter(sfari_raw, is.na(ensembl_id)) %>%
    left_join(missing_id_matches, by = c("gene_symbol" = "gene_name")) %>%
    dplyr::select(-ensembl_id) %>%
    dplyr::select(status, gene_symbol, gene_name, "ensembl_id" = gene_id, everything()) %>%
    filter(!is.na(ensembl_id))
  ) 

sfari_cleaned %>% nrow() # 1092 human genes

## Convert to Dmel genes
fb_diopt <- read_tsv("./data/Flybase/dmel_human_orthologs_disease_fb_2019_03.tsv", skip = 4) %>%
  janitor::clean_names() %>%
  left_join(hsap_tx2gene, by = c("human_gene_symbol" = "gene_name")) 

hsap_dmel_conversion <- fb_diopt %>%
  filter(diopt_score >= 8) %>%
  dplyr::select("dmel_gene_id" = number_number_dmel_gene_id, dmel_gene_symbol, human_gene_hgnc_id, human_gene_symbol, diopt_score, "human_gene_id" = gene_id)

sfari_dmel <- sfari_cleaned %>%
  left_join(hsap_dmel_conversion, by = c("ensembl_id" = "human_gene_id"))

sfari_dmel %>% write_csv("./data/SFARI/sfari_left_join_dmel_diopt.csv")

sfari_dmel_genes <- sfari_dmel %>%
  pull(dmel_gene_id) %>% unique()

length(sfari_dmel_genes) # 632 dmel genes converted from sfari

## Compare with 1700 
source("./Rscript/run_hypergeometric_test.R")

intersection_group <- intersect(sfari_dmel_genes, id_interest$dmel_gene_id) %>% length() # 229 genes
glial_group <- nrow(id_interest)
sfari_group <- intersect(sfari_dmel_genes, gl_nd$dmel_gene_id) %>% length()
total_group <- nrow(gl_nd)

run_hypergeometric_test(intersection_group, glial_group, sfari_group, total_group)






View(sfari_dmel)






View(id_interest)




