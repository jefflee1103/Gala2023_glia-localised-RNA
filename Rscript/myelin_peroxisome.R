# ------ For Ilan/Julia
## Find myelin-localised trasncripts involved in peroxisome
## Re-analysis of Thakurela 2016 (Mmus)

# ------ Create TPM dataframe
library(tidyverse)
library(qs)

Mmus_avgTPM_tidy <-readRDS("./data/RNAseq/quant_results/Mmus_kallisto_TPM_tidy.RDS")

myelin_avgTPM_tidy <- Mmus_avgTPM_tidy %>%
  filter(str_detect(library_group, "Myelin")) %>%
  group_by(library_group, gene_id) %>%
  summarise(avg_TPM = mean(TPM))

myelin_avgTPM_wide <- myelin_avgTPM_tidy %>%
  pivot_wider(values_from = avg_TPM, names_from = library_group)

myelin_present_filtered <- myelin_avgTPM_wide %>%
  filter(if_any(contains("txn"), ~ .x > 10))

# ------ Add BiomaRt GO data
library(biomaRt)

mouse_mart96 <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", 
                        host = "https://apr2019.archive.ensembl.org")
listAttributes(mouse_mart96) -> mouse_biomart_attributes
id_dictionary <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "go_id", "name_1006", "namespace_1003", "hsapiens_homolog_associated_gene_name"),
                       mart = mouse_mart96,
                       filters = "ensembl_gene_id",
                       values = myelin_present_filtered$gene_id)

human_orthologs <- getBM(attributes = c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"),
                         mart = mouse_mart96,
                         filters = "ensembl_gene_id",
                         values = myelin_present_filtered$gene_id)
human_orthologs <- human_orthologs %>%
  group_by(ensembl_gene_id) %>%
  summarise(human_homologs = str_c(hsapiens_homolog_associated_gene_name, collapse = " | "))

genes_involved_in_peroxisome_function <- id_dictionary %>%
  filter(str_detect(name_1006, "peroxisome")) %>%
  pull(ensembl_gene_id) %>%
  unique()

gene_names <- id_dictionary %>%
  dplyr::select(ensembl_gene_id, external_gene_name, description) %>%
  distinct()

myelin_present_filtered_anno <- myelin_present_filtered %>%
  left_join(gene_names, by = c("gene_id" = "ensembl_gene_id")) %>%
  left_join(human_orthologs, by = c("gene_id" = "ensembl_gene_id")) %>%
  mutate(GO_peroxisome_function = if_else(gene_id %in% genes_involved_in_peroxisome_function, "Yes", "No")) %>%
  dplyr::select("ensembl_gene_id" = gene_id, "ensembl_gene_name" = external_gene_name, human_homologs, GO_peroxisome_function, "full_name" = description, everything()) %>%
  filter(!str_detect(full_name, "predicted gene")) %>%
  arrange(desc(GO_peroxisome_function), desc(`Myelin-24mo_protrusion_txn`))

write_csv(myelin_present_filtered_anno, "./output/analysis/peroxisome/myelin_localised_transcripts_df.csv")

# ----- Summary of peroxisome GO genes
go_summary_df <- id_dictionary %>%
  filter(str_detect(name_1006, "peroxisome")) %>%
  left_join(human_orthologs) %>%
  group_by(name_1006, namespace_1003) %>%
  summarise(
    gene_count = n(), 
    mouse_genes = str_c(external_gene_name, collapse = " | "),
    human_genes = str_c(human_homologs, collapse = " | ")
    ) %>%
  arrange(desc(gene_count)) %>%
  dplyr::rename("GO_entry" = name_1006, "GO_group" = namespace_1003)

write_csv(go_summary_df, "./output/analysis/peroxisome/myelin_localised_transcripts_peroxisome_go.csv")
