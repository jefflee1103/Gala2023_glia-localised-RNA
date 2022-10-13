# Identify if there are any enrichment of human disease assosication 
# from glia-localised transcript data using flybase annotations 

# ----- Environment 

library(tidyverse)
library(ggrepel)
library(qs)
library(furrr)
plan(multisession, workers = 8)

# ----- Prepare Flybase disease association dataframe 

## FB disease model annotations
dismodel_raw <- read_tsv("./data/Flybase/disease_model_annotations_fb_2022_04.tsv", skip = 4) %>%
  setNames(c("current_gene_id", "gene_symbol", "hgnc_id", "do_qualifier", "do_id", "do_term", "allele_used_in_model_fbal", "allele_used_in_model_symbol", "orthology_with_hgnc_id", "orthology_with_symbol", "evidence", "reference")) %>%
  filter(!str_detect(gene_symbol, "[:alnum:]+\\\\")) %>%
  filter(!str_detect(do_qualifier, "DOES NOT"))

dismodel_cleaned <- dismodel_raw %>%
  mutate(do_term = str_replace(do_term, " type [:alnum:]+$", "")) %>%
  mutate(do_term = str_replace(do_term, " group [:alnum:]+$", "")) %>%
  mutate(do_term = str_replace(do_term, " variant type$", "")) %>%
  mutate(do_term = str_replace(do_term, " [:digit:]+$", "")) %>%
  mutate(do_term = str_replace(do_term, " [:digit:]+[:alpha:]$", "")) %>%
  mutate(do_term = str_replace(do_term, "-[:digit:]+$", "")) %>%
  mutate(do_term = str_replace(do_term, " [IVX]+$", "")) %>%
  mutate(do_term = str_replace(do_term, " [IVX]+[:alnum:]$", "")) %>%
  mutate(do_term = str_replace(do_term, "hypogonadotropic hypogonadism .+", "hypogonadotropic hypogonadism")) %>%
  mutate(do_term = str_replace(do_term, "mitochondrial complex .+", "mitochondrial complex deficiency")) %>%
  mutate(do_term = str_replace(do_term, "syndromic X-linked intellectual .+", "syndromic X-linked intellectual disability")) %>%
  mutate(do_term = str_replace(do_term, "syndromic X-linked mental retardation .+", "syndromic X-linked mental retardation")) %>%
  mutate(do_term = str_replace(do_term, "short-rib thoracic dysplasia .+", "short-rib thoracic dysplasia")) %>%
  mutate(do_term = str_replace(do_term, "neurodevelopmental disorder .+", "neurodevelopmental disorder")) %>%
  mutate(do_term = str_replace(do_term, "mucolipidosis .+", "mucolipidosis")) %>%
  mutate(do_term = str_replace(do_term, "microcephaly .+", "microcephaly")) %>%
  mutate(do_term = str_replace(do_term, "lissencephaly .+", "lissencephaly")) %>%
  mutate(do_term = str_replace(do_term, "Leber congenital amaurosis .+", "Leber congenital amaurosis")) %>%
  mutate(do_term = str_replace(do_term, "immunodeficiency .+", "immunodeficiency")) %>%
  mutate(do_term = str_replace(do_term, "spinal muscular atrophy .+", "spinal muscular atrophy")) %>%
  mutate(do_term = str_replace(do_term, "hypomyelinating leukodystrophyy .+", "hypomyelinating leukodystrophyy")) %>%
  mutate(do_term = str_replace(do_term, "spinal muscular atrophy .+", "spinal muscular atrophy"))

## Resolve current gene_id to ens99 gene_id
ens99 <- read_csv("./data/Flybase/Dmel_tx2gene_ENSEMBL_v99.csv") %>%
  dplyr::select(gene_id, gene_name) %>% distinct() %>%
  mutate(ens99_gene_id = gene_id)
syno <- read_tsv("./data/Flybase/fbgn_annotation_ID_fb_2022_04.tsv", skip = 4) %>%
  setNames(c("gene_symbol", "organism", "current_gene_id", "secondary_gene_ids", "annotation_id", "secondary_annotation_ids")) %>%
  filter(organism == "Dmel") %>%
  dplyr::select(current_gene_id, secondary_gene_ids)

resolve_tmp <- left_join(dismodel_cleaned, ens99, by = c("current_gene_id" = "gene_id"))
resolve_noconflict <- filter(resolve_tmp, !is.na(ens99_gene_id)) %>%
  dplyr::select(colnames(dismodel_cleaned), "gene_id" = current_gene_id)
resolve_yesconflict <- filter(resolve_tmp, is.na(ens99_gene_id)) %>%
  left_join(syno, by = "current_gene_id") %>%
  filter(!is.na(secondary_gene_ids)) %>%
  separate_rows(secondary_gene_ids, sep = ",") %>% 
  distinct() %>%
  mutate(is_present_ens99 = if_else(secondary_gene_ids %in% ens99$gene_id, TRUE, FALSE))

resolve_yesconflict_verdict <- resolve_yesconflict %>%
  filter(is_present_ens99 == TRUE) %>%
  dplyr::select(colnames(dismodel_cleaned), -current_gene_id, "gene_id" = secondary_gene_ids) %>%
  dplyr::select(gene_id, everything())

dismodel_ens99 <- bind_rows(resolve_noconflict, resolve_yesconflict_verdict) %>%
  left_join(ens99, by = "gene_id") %>%
  dplyr::select(-gene_symbol, -ens99_gene_id) %>%
  dplyr::select(gene_id, gene_name, everything())

# write_csv(dismodel_ens99, "./data/Flybase/disease_model_annotations_fb_2022_04_ens99.csv")

# ----- Perform statistics of gene set enrichment

## Prepare dataframe for statstics 
gl_nd <- qread("./output/glia-localised-neurodegenerative-disease-incorporated.qs")
id_interest <- gl_nd %>%
  filter(rna_in_protrusion >= 8) %>%
  filter(`subperineurial glial cell` == TRUE | 
           `perineurial glial sheath` == TRUE |
           `adult brain perineurial glial cell` == TRUE |
           `ensheathing glial cell` == TRUE)

gl_disease_count <- unique(dismodel_ens99$do_term) %>%
  set_names() %>%
  future_map_dfr(~{
    do_genes <- filter(dismodel_ens99, do_term == .x) %>% pull(gene_id) %>% unique()
    tibble(
      count_in_genome = intersect(ens99$gene_id, do_genes) %>% length(),
      count_in_gl = intersect(id_interest$dmel_gene_id, do_genes) %>% length()
    )
  }, .id = "disease")

total_background_count <- length(ens99$gene_id)
total_target_count <- length(id_interest$dmel_gene_id)

gl_disease_statistics <- gl_disease_count %>%
  mutate(
    pvalue = phyper(
      count_in_gl - 1, 
      count_in_genome, 
      total_background_count - count_in_genome,
      total_target_count,
      lower.tail = FALSE
      )
    ) %>%
  mutate(expected_count = total_target_count * (count_in_genome / total_background_count)) %>%
  mutate(foldchange = count_in_gl / expected_count) %>%
  mutate(log2foldchange = log2(foldchange)) %>%
  mutate(adj_pvalue = p.adjust(pvalue, method = "bonferroni"))

## Dataframe for list of genes 
gl_disease_gene_list <- gl_disease_statistics %>%
  dplyr::select(disease, log2foldchange, adj_pvalue, count_in_genome,"count_in_localised_genes" = count_in_gl) %>%
  filter(-log(adj_pvalue) > 10) %>%
  arrange(adj_pvalue) %>%
  mutate(gene_list = map_chr(disease, ~ {
    filter(dismodel_ens99, do_term == .x) %>% 
      filter(gene_id %in% id_interest$dmel_gene_id) %>%
      pull(gene_name) %>% unique() %>% paste(collapse = ", ") %>% sort()
  })) 

write_tsv(gl_disease_gene_list, "./output/analysis/gl_disease_gene_list_top-pvalue.txt")

# ----- Plot

gl_disease_statistics %>%
  filter(log2foldchange != -Inf) %>%
  mutate(put_label = if_else(adj_pvalue < 0.01, disease, NA_character_)) %>%
  ggplot(aes(x = log2foldchange, y = -log(adj_pvalue), colour = count_in_gl)) +
  geom_point(alpha = 1, stroke = 0, size = 3) +
  geom_text_repel(aes(label = put_label), 
    max.overlaps = 100, hjust = 0,
    cex = 2.2, colour = "gray40") +
  geom_hline(yintercept = -log(0.01), linetype = "dashed", colour = "gray80") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray80") +
  labs(
    title = "Enrichment of disease-associated genes within glial protrusion-localised transcripts",
    subtitle = "Vocalno plot: Fisher's exact test against full fly genome",
    x = "log2FoldChange",
    y = "Bonferroni-adjusted p-value (-log10)",
    colour = "Localised genes"
  ) + 
  scale_colour_viridis_c() +
  coord_cartesian(xlim = c(-3.5, 3.5)) +
  theme_classic(base_size = 10) +
  annotate(geom = "text", x = -3, y = -log(0.01) + 2, label = "p=0.01", cex = 3, alpha = 0.5)

ggsave("./output/graphics/disease_association_volcano_plot.pdf", 
  width = 10, height = 8, useDingbats = FALSE)


















for(i in (conflict$dmel_gene_name)) {
  id_interest <- id_interest %>%
    mutate(dmel_gene_id = 
      if_else(dmel_gene_name == i,
        filter(ens99, gene_name == i) %>% pull(ens99_gene_id),
        dmel_gene_id))
}


conflict <- id_interest %>%
  left_join(ens99, by = c("dmel_gene_id" = "gene_id")) %>%
  filter(is.na(ens99_gene_id))

View(ens99)

syno


intersect(id_interest$dmel_gene_id, ens99$gene_id) %>% length()

id_interest %>% write_tsv("~/Desktop/glia-protrusion-localised-id-interest.txt")
