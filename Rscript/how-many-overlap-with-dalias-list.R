library(tidyverse)
library(qs)

## Get gl + genome data
ens99 <- read_csv("./data/Flybase/Dmel_tx2gene_ENSEMBL_v99.csv") %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct()
gl_nd <- qread("./output/glia-localised-neurodegenerative-disease-incorporated.qs")
id_interest <- gl_nd %>%
  filter(rna_in_protrusion >= 8) %>%
  filter(`subperineurial glial cell` == TRUE | 
           `perineurial glial sheath` == TRUE |
           `adult brain perineurial glial cell` == TRUE |
           `ensheathing glial cell` == TRUE)

## Get Dalia's list data 
dalia_list <- tibble(
  X1 = c(
    "nrv2",
    "Cip4",
    "Lac",
    "Nrg",
    "Gs2",
    "cold",
    "Atpalpha",
    "lost",
    "CG1648",
    "CG42342",
    "kst",
    "Pdi",
    "ORMDL",
    "sdk",
    "alpha-Cat",
    "Vha55",
    "Flo2",
    "shot",
    "Nrx-IV",
    "Rab11"
  )
)

dalia_list_v2 <- tibble(
  X1 = c(
    "alpha-Cat",
    "arm",
    "Atpalpha",
    "CG1648",
    "CG42342",
    "Cip4",
    "cno",
    "cold",
    "Flo2",
    "Gli",
    "lost",
    "kst",
    "l(1)G0320",
    "Lac",
    "Nrg",
    "nrv2",
    "Nrx-IV",
    "Pdi",
    "Sdc",
    "Vha55",
    "zip"
  )
)

gene_list <- dalia_list_v2 %>%
  left_join(ens99, by = c("X1" = "gene_name")) %>%
  dplyr::select(gene_id, "gene_name" = X1) %>%
  mutate(
    is_mouse_dmel_convertible =
      if_else(gene_id %in% gl_nd$dmel_gene_id, TRUE, FALSE)
  ) %>%
  mutate(
    is_in_1700 =
      if_else(gene_id %in% id_interest$dmel_gene_id, TRUE, FALSE)
  ) %>%
  left_join(id_interest, by = c("gene_id" = "dmel_gene_id")) %>%
  arrange(desc(is_mouse_dmel_convertible), desc(is_in_1700))

View(gene_list)
























"ORMDL" %in% id_interest$dmel_gene_name

"FBgn0013997" %in% id_interest$dmel_gene_id

View(id_interest)
View(gl_nd)
