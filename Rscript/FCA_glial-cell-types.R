library(SCopeLoomR)
library(tidyverse)
library(qs)
library(furrr)
library(patchwork)

plan(multisession, workers = 4)

## loom
loom_path <- '~/Downloads/s_fca_biohub_glial_cell_10x.loom'
loom <- open_loom(loom_path)

embeddings <- get_embeddings(loom) %>%
  imap_dfc(~ pluck(.x) %>% 
             as.data.frame() %>%
             setNames(c(paste0(.y, "_x"), paste0(.y, "_y")))) %>%
  rownames_to_column(var = "cell_id") %>%
  janitor::clean_names()

cell_anno <- get_cell_annotation(loom) %>%
  rownames_to_column(var = "cell_id") %>%
  left_join(embeddings, by = "cell_id") %>%
  dplyr::select(cell_id, age, annotation, annotation__ontology_id, n_counts, n_genes, percent_mito, sex, tissue, contains("hvg"))

cell_anno_summary <- cell_anno %>%
  dplyr::select(annotation, annotation__ontology_id) %>%
  group_by(annotation, annotation__ontology_id) %>%
  summarise(n_cell = n()) 

## Glia annotations to keep
glia_anno <- c(
  "adult antenna glial cell",
  "adult brain cell body glial cell",
  "adult brain perineurial glial cell",
  "adult glial cell",
  "adult lamina epithelial/marginal glial cell",
  "adult optic chiasma glial cell",
  "adult reticular neuropil associated glial cell",
  "cell body glial cell",
  "CNS surface associated glial cell",
  "ensheathing glial cell",
  "optic-lobe-associated cortex glial cell",
  "perineurial glial sheath",
  "peripheral glial cell",
  "subperineurial glial cell"
)

glia_easy_anno <- tibble(
  annotation = glia_anno,
  easy_annotation = c(
    "antenna glia",
    "uncategorised 01",
    "perineurial glia",
    "uncategorised 02",
    "marginal glia",
    "optic chiasma glia",
    "astrocyte-like glia",
    "uncategorised 03",
    "uncategorised 04",
    "ensheathing glia",
    "cortex glia",
    "perineurial glia",
    "uncategorised 05 - peripheral",
    "subperineurial glia"
  )
  )
  

## Filter and get pure-glia data
all_glia_cell_anno <- filter(cell_anno, annotation %in% glia_anno) %>%
  left_join(glia_easy_anno)
all_glia_cell_ids <- pull(all_glia_cell_anno, cell_id)
# all_glia_dgem <- filter(dgem_tidy, cell_id %in% all_glia_cell_ids)
# qsave(all_glia_dgem, "~/Downloads/all_glia_dgem.qs")

rm(dgem, dgem_tidy, embeddings, cell_anno, cell_anno_summary)

# all_glia_dgem_by_gene <- all_glia_dgem %>%
#   group_by(gene_name) %>%
#   summarise(n_exp_cell = sum(dgem > 0),
#             sum_dgem = sum(dgem)) %>%
#   ungroup() %>%
#   mutate(pct_exp_cell = (n_exp_cell / length(unique(all_glia_dgem$cell_id))) * 100) %>%
#   mutate(exp_level_within_expressing_cell = sum_dgem / n_exp_cell) 

## Glia dgem
all_glia_dgem <- qread("~/Downloads/all_glia_dgem.qs")

all_glia_dgem_by_gene <- qread("./data/FCA/glia_type_dgem_by_gene.qs") %>%
  left_join(glia_easy_anno, by = c("glia_cell_type" = "annotation"))

glia_type_dgem_by_gene_easy <- glia_easy_anno$easy_annotation %>%
  set_names() %>%
  map_dfr(function(x){
    cell_type_id <- filter(all_glia_cell_anno, easy_annotation == x) %>% pull(cell_id) %>% unique()
    filter(all_glia_dgem, cell_id %in% cell_type_id) %>%
      group_by(gene_name) %>%
      summarise(n_exp_cell = sum(dgem > 0),
                sum_dgem = sum(dgem)) %>%
      ungroup() %>%
      mutate(pct_exp_cell = (n_exp_cell / length(cell_type_id)) * 100) %>%
      mutate(exp_level_within_expressing_cell = sum_dgem / n_exp_cell) %>%
      mutate(glia_cell_type_easy = x)
  }, .progress = TRUE)

# glia_type_dgem_by_gene_easy_withtotal <- glia_type_dgem_by_gene_easy %>%
#   bind_rows(
#     (glia_type_dgem_by_gene_easy %>%
#        mutate(glia_cell_type_easy = "All subtypes"))
#   )


## Gene counts per glia type
all_glia_cell_anno %>%
  ggplot(aes(x = easy_annotation, y = n_genes, colour = easy_annotation)) +
  ggbeeswarm::geom_quasirandom(size = 1, alpha = 0.3, stroke = 0) +
  scale_colour_viridis_d() +
  labs(title = "Number of genes detected per cell",
       subtitle = "Fly cell atlas - single-nucleus RNA-seq",
       x = "",
       y = "Gene count") + 
  theme_gray(base_size = 12) +
  theme(legend.position = "none") +
  coord_flip() +
  ggbeeswarm::geom_quasirandom(aes(x = "TOTAL", y = n_genes),
                               size = 1, alpha = 0.3)



## summary 

gl_dmel <- qread("./output/glia-localised-rna_dmel-converted.qs")
gl_fca <- qread("./output/glia-localised-rna_fca-incorporated.qs")
present_genes <- gl_dmel %>%
  filter(rna_in_protrusion >= 8) %>%
  pull(dmel_gene_name)
id_interest_1700_genes <- gl_fca %>%
  filter(rna_in_protrusion >= 8) %>%
  filter(`subperineurial glial cell` == TRUE | 
           `perineurial glial sheath` == TRUE |
           `adult brain perineurial glial cell` == TRUE |
           `ensheathing glial cell` == TRUE) %>% pull(dmel_gene_name)


gliatype_summary_df <- glia_type_dgem_by_gene_easy %>%
  distinct() %>%
  mutate(is_RNA_present = if_else(gene_name %in% present_genes, TRUE, FALSE)) %>% 
  mutate(is_id_1700 = if_else(gene_name %in% id_interest_1700_genes, TRUE, FALSE)) %>%
  group_by(glia_cell_type_easy) %>%
  summarise(
    scRNAseq_expressed = sum(pct_exp_cell > 2.5),
    RNA_present = sum(pct_exp_cell > 2.5 & is_RNA_present == TRUE),
    overlap_1700 = sum(pct_exp_cell > 2.5 & is_id_1700 == TRUE)
  )

all_types_scRNAseq_expressed <- glia_type_dgem_by_gene_easy %>%
    filter(pct_exp_cell > 2.5) %>% pull(gene_name) %>% unique()

allglia_summary_df <- tibble(
  glia_cell_type_easy = "All Glia",
  scRNAseq_expressed = length(all_types_scRNAseq_expressed),
  RNA_present = intersect(all_types_scRNAseq_expressed, present_genes) %>% length(),
  overlap_1700 = 1700
)

summary_df_long <- bind_rows(allglia_sumamry_df, gliatype_summary_df) %>%
  pivot_longer(cols = -glia_cell_type_easy, names_to = "type", values_to = "count") %>%
  mutate(type = case_when(
    type == "scRNAseq_expressed" ~ "1. scRNAseq expressed",
    type == "RNA_present" ~ "2. Localised RNA present (8/12 libraries)",
    type == "overlap_1700" ~ "3. Overlap with 1700",
  ))

summary_df_long %>%
  mutate(type = fct_rev(type)) %>%
  mutate(glia_cell_type_easy = fct_rev(glia_cell_type_easy)) %>%
  ggplot(aes(x = glia_cell_type_easy, y = count, fill = type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(x = glia_cell_type_easy, y = count + 200, label = scales::comma(count), group = type),
             position = position_dodge(width = 0.8), cex = 3.5, inherit.aes = FALSE) + 
  labs(title = "Expressed and localised genes per glial subtype",
       subtitle = "From Fly Cell Atlas, Astrocytes = alrm-expressing cells", 
       x = "",
       y = "Gene count",
       fill = "") + 
  scale_y_continuous(label = scales::comma) + 
  scale_fill_manual(values = rev(c("gray70", "goldenrod3", "goldenrod1"))) + 
  theme_minimal(base_size = 15) + 
  theme(legend.position = "bottom") + 
  coord_flip() +
  guides(
    fill = guide_legend(direction = "vertical", reverse = TRUE)
  )

ggsave("~/Desktop/expressed_localised_genes_per_glial_subtype.jpg", height = 11, width = 9)


## overlap with Maria's AI data
maria_ai <- read_csv("~/Desktop/maria_ai.csv")
fca_glia_genenames <- tibble(gene_names = unique(glia_type_dgem_by_gene_easy$gene_name))

maria_ai %>%
  left_join(glia_type_dgem_by_gene_easy, by = c("maria_ai" = "gene_name")) %>%
  mutate(expressing_glia = if_else(pct_exp_cell > 2.5, glia_cell_type_easy, NA)) %>%
  group_by(maria_ai) %>%
  summarise(
    expressing_glia = paste(sort(unique(expressing_glia)), collapse = "; ")
  ) %>%
  write_tsv("~/Desktop/maria_ai66_expressing_glia.txt")

test
maria_ai$maria_ai %>% unique() %>% length()

test$maria_ai %>% unique() %>% length()

setdiff(maria_ai$maria_ai, test$maria_ai)





##
plot_gene_tsne <- function(gene, palette, limit_counts){
  counts <- all_glia_dgem %>% filter(gene_name == gene)
  df <- left_join(all_glia_cell_anno, counts, by = "cell_id") %>%
    mutate(count = ifelse(dgem == 0, NA, dgem))
  ggplot(df, aes(x = hvg_t_sne_x, y = hvg_t_sne_y, colour = count)) +
    geom_point(size = 1, alpha = 0.9, stroke = 0) +
    labs(title = paste0(gene, " expression"),
         subtitle = "Fly cell atlas - single-nucleus RNA-seq") + 
    # scale_colour_gradientn(colours = c("gray90", "gray90", palette),
    #                        values = c(0, 0.1, 1),
    #                        limits = c(0, limit_counts)) +
    scale_colour_gradient(low = "gray50", high = palette, na.value = "gray80",
                          limits = c(0, limit_counts),
                          oob = scales::squish) +
    theme_void(base_size = 15) +
    theme(legend.position = "bottom")
}

plot_gene_tsne("nrv2", "red", 3)

















