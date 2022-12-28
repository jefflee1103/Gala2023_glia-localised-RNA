# Identify if there are any enrichment of human disease assosication 
# from glia-localised transcript data using flybase annotations 

# ----- Environment 

library(tidyverse)
library(ggrepel)
library(qs)
library(furrr)
plan(multisession, workers = 4)

# ----- Prepare Flybase disease association dataframe 

## FB disease model annotations
dismodel_raw <- read_csv("./data/REACTOME/flymine_gene_pathway_level1.csv") %>% 
  setNames(c("pathways_name_level0", "gene_secondary_identifier", "gene_symbol", "pathways_identifier_level1", "pathways_name_level1", "datasets_name","current_gene_id")) 

#%>% filter(!str_detect(do_qualifier, "DOES NOT"))

dismodel_cleaned <- dismodel_raw 

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


dismodel_ens99 %>%
  group_by(pathways_name_level1) %>%
  summarise(count = n()) -> tmp


# ----- Perform statistics of gene set enrichment

## Prepare dataframe for statistics 
gl_nd <- qread("./output/glia-localised-neurodegenerative-disease-incorporated.qs")
id_interest <- gl_nd %>%
 filter(rna_in_protrusion >= 8) %>%
 filter(`subperineurial glial cell` == TRUE |
          `perineurial glial sheath` == TRUE |
          `adult brain perineurial glial cell` == TRUE |
          `ensheathing glial cell` == TRUE)

gl_disease_count <- unique(dismodel_ens99$pathways_name_level1) %>%
  set_names() %>%
  future_map_dfr(~{
    do_genes <- filter(dismodel_ens99, pathways_name_level1 == .x) %>% pull(gene_id) %>% unique()
    tibble(
      count_in_genome = intersect(ens99$gene_id, do_genes) %>% length(),
      count_in_gl = intersect(id_interest$dmel_gene_id, do_genes) %>% length()
    )
  }, .id = "name")

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
  dplyr::select(name, log2foldchange, adj_pvalue, count_in_genome,"count_in_localised_genes" = count_in_gl) %>%
  filter(-log(adj_pvalue) > 10) %>%
  arrange(adj_pvalue) %>%
  mutate(gene_list = map_chr(name, ~ {
    filter(dismodel_ens99, pathways_name_level1 == .x) %>% 
      filter(gene_id %in% id_interest$dmel_gene_id) %>%
      pull(gene_name) %>% unique() %>% paste(collapse = ", ") %>% sort()
  })) 

#write_tsv(gl_disease_gene_list, "./output/analysis/gl_")

# ----- Plot

plottingdf <- gl_disease_statistics %>%
  filter(log2foldchange != -Inf) %>%
  mutate(put_label = if_else(-log10(adj_pvalue) > 15 | log2foldchange > 2.5, name, NA_character_)) %>%
  mutate(put_label = if_else(-log10(adj_pvalue) < 5, NA_character_, put_label))
  
plottingdf$name2 <- dismodel_ens99$pathways_name_level0[match(plottingdf$name, dismodel_ens99$pathways_name_level1)]

plottingdf <- plottingdf %>%
  group_by(name2) %>%
  mutate(level1_occurrence = sum(adj_pvalue < 0.01)) %>%
  ungroup() %>%
  mutate(name2 = if_else(level1_occurrence < 2, "other", name2))

level1_order <- plottingdf %>%
  arrange(desc(level1_occurrence)) %>%
  pull(name2) %>%
  unique()

colours <- c(colorspace::qualitative_hcl(palette = "Dynamic", n = 11), "gray70")

plottingdf %>%
  mutate(name2 = fct_relevel(name2, level1_order)) %>%
  ggplot(aes(x = log2foldchange, y = -log10(adj_pvalue), colour = name2, size = count_in_gl)) +
  geom_point(alpha = 0.9, stroke = 0) +
  geom_text_repel(aes(label = put_label), 
    max.overlaps = 100, hjust = 0,
    cex = 1, colour = "gray40", ) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", colour = "gray80", size = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray80", size = 0.3) +
  labs(
    x = "log2FoldChange",
    y = "Bonferroni-adjusted p-value (-log10)",
    colour = "Reactome group (n>2 significant)",
    size = "Localised genes"
  ) + 
  scale_colour_manual(values = colours) +
  coord_cartesian(xlim = c(0, 3.5)) +
  theme_classic(base_size = 6) +
  annotate(geom = "text", x = 0.25, y = -log(0.01) + 2, label = "p=0.01", cex = 1.5, alpha = 0.5) +
  guides(colour = guide_legend(keywidth = 0.3, keyheight = 0.5),
  size = guide_legend(keywidth = 0.3, keyheight = 0.5))

ggsave("./output/graphics/reactome_pathway_association_volcano_plot.pdf",
  width = 13.5 * 0.3937, height = 10 * 0.3937,
  device = cairo_pdf
)



