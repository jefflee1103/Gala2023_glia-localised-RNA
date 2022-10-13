# ----- Environment

library(tidyverse)
library(readxl)
library(qs)

# ----- Import data and clean up

## Get screen paper data
cpti_figure_id <- read_csv("./data/screen-paper/screen_cpti-figureid-mapper.csv")
cpti_figure_id_tidy <- cpti_figure_id %>%
  dplyr::select(-overview_id) %>%
  dplyr::rename_with(~str_replace(.x, "_[:lower:]+", "_figureid"), contains("_id", ignore.case = FALSE)) %>%
  dplyr::rename("gene_id" = FBgn_ID) %>%
  mutate(across(contains("figureid"), ~ as.character(.x))) %>%
  pivot_longer(cols = contains("_figureid"), names_to = "figure_type", values_to = "figure_id")


### check if cpti fbgns are compatible with ensembl 99 
ens99 <- read_csv("./data/Flybase/Dmel_tx2gene_ENSEMBL_v99.csv") %>%
  dplyr::select(gene_id, gene_name) %>% distinct()

intersect(unique(cpti_figure_id$FBgn_ID), ens99$gene_id) %>% length()


## Get Josh's glial process screen data 
nmj_glia <- read_excel("./data/screen-paper/screen-paper_supp3.xlsx", sheet = 1)
cb_glia <- read_excel("./data/screen-paper/screen-paper_supp3.xlsx", sheet = 4)
mb_glia <- read_excel("./data/screen-paper/screen-paper_supp3.xlsx", sheet = 5)

glia_rna_present <- list(nmj_glia, cb_glia, mb_glia) %>%
  set_names(c("nmj", "cb", "mb")) %>%
  map(~ janitor::clean_names(.x) %>% 
        dplyr::select(figure_id, contains("rna_in")) %>%
        filter(if_any(contains("rna_in"), function(x){x == "yes"})))

glia_rna_present %>%
  map(~ nrow(.x))

glia_rna_present_genes <- glia_rna_present %>%
  map(~ mutate(.x, figure_id = as.character(figure_id)) %>%
        left_join(cpti_figure_id_tidy, by = "figure_id") %>% 
        pull(gene_id)) %>%
  reduce(append) %>%
  unique()

length(glia_rna_present_genes)



