library(tidyverse)
library(qs)

gl_nd <- qread("./output/glia-localised-neurodegenerative-disease-incorporated.qs")
id_interest <- gl_nd %>%
  filter(rna_in_protrusion >= 8) %>%
  filter(`subperineurial glial cell` == TRUE | 
           `perineurial glial sheath` == TRUE |
           `adult brain perineurial glial cell` == TRUE |
           `ensheathing glial cell` == TRUE)


dalia_list <- read_delim(pipe("pbpaste"), delim = "/t", col_names = FALSE)

ens99 <- read_csv("./data/Flybase/Dmel_tx2gene_ENSEMBL_v99.csv") %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct()

gene_list <- left_join(dalia_list, ens99, by = c("X1" = "gene_name")) %>%
  dplyr::select(gene_id, "gene_name" = X1) %>%
  mutate(is_in_1700 = if_else(gene_id %in% id_interest$dmel_gene_id, TRUE, FALSE)) %>%
  arrange(desc(is_in_1700), gene_name)


gene_list







dalia_list
head(id_interest)








"ORMDL" %in% id_interest$dmel_gene_name

"FBgn0013997" %in% id_interest$dmel_gene_id

View(id_interest)
