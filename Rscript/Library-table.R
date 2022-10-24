
# Create glial protrusion library table -------------------------------------------------------

library(tidyverse)
library(readxl)
library(qs)
library(rcartocolor)

# Get plotting dataframe ----------------------------------------------------------------------

## Glial data
avgTPM_wide <- readRDS("./data/RNAseq/quant_results/all-glia_avgTPM_wide.RDS")

### Only keep the 12 libraries that we are interested in 
protrusion_RNA_library_info <- read_csv("./data/RNAseq/protrusion_RNA_library_info.csv") %>%
  unite(col = "index", c(Model.system, Species, Separation.method, Data.type, Reference), sep = "\n\n")
protrusion_RNA_libraries <- protrusion_RNA_library_info$library
avgTPM <- avgTPM_wide %>%
  dplyr::select(gene_id, protrusion_RNA_libraries) %>%
  pivot_longer(cols = contains(c("trap", "txn")),
               names_to = "library",
               values_to = "TPM")

## Neurite data
neurite_enriched <- read_excel("./data/RNAseq/external_data/Supplementary online tables_neurite-enriched.xlsx",
                               skip = 1) %>% janitor::clean_names()

neurite_present <- read_excel("./data/RNAseq/external_data/Supplementary online tables_Chekulaeva_most-abundant.xlsx",
                              skip = 1) %>% janitor::clean_names()

## Annotate and summarise
### Genes that are commonly detected in at least 3 datasets
commonly_expressed_genes <- avgTPM %>%
  group_by(gene_id) %>%
  summarise(n_exp_lib = sum(TPM > 10)) %>% ungroup() %>%
  filter(n_exp_lib >= 3) %>%
  pull(gene_id)

avgTPM_anno <- avgTPM %>%
  mutate(protrusion_RNA = if_else(TPM > 10, TRUE, FALSE)) %>%
  mutate(is_common = if_else(gene_id %in% commonly_expressed_genes, "Detected in ≥ 3 datasets", "Detected in < 3 datasets")) %>%
  mutate(is_neurite_present = if_else(gene_id %in% neurite_present$gene_id, "Shared with neurites", "Glia only"))

avgTPM_summary_glia_dataset <- avgTPM_anno %>%
  filter(protrusion_RNA == TRUE) %>%
  group_by(library, is_common) %>%
  summarise(count = n()) %>% ungroup() %>%
  left_join(protrusion_RNA_library_info)

avgTPM_summary_glia_neurite <- avgTPM_anno %>%
  filter(protrusion_RNA == TRUE) %>%
  group_by(library, is_neurite_present) %>%
  summarise(count = n()) %>% ungroup() %>%
  left_join(protrusion_RNA_library_info)

expression("">=3)

# Plot overall study summary ------------------------------------------------------------------

colour <- carto_pal(12, "Vivid")

c("#E58606", "#5D69B1", "#52BCA3", "#99C945", "#CC61B0", "#24796C", "#DAA51B", "#2F8AC4", "#764E9F", "#ED645A", "#CC3A8E", "#A5AA99")

colour_to_use <- c("gray75", "#E58606")

avgTPM_summary_glia_dataset %>%
  mutate(index = str_replace(index, "trans", "\ntrans")) %>%
  mutate(index = str_replace(index, "al,", "al,\n")) %>%
  # mutate(is_common, str_replace(is_common, "≥", "\u2265")) %>%
  ggplot(aes(
    x = index,
    y = count, 
    fill = is_common
  )) +
  labs(fill = "",
       x = "",
       y = "Detected gene count (TPM >10)") + 
  geom_col(width = 0.6, alpha = 0.9) + 
  scale_fill_manual(values = colour_to_use) +
  scale_y_continuous(expand = c(0, 100)) + 
  theme_classic(base_size = 6) + 
  theme(legend.position = "right",
    axis.text.x = element_text(size = 3),
    legend.text = element_text(size = 4)) +
  guides(fill = guide_legend(
    reverse = FALSE,
    keywidth = 0.3,
    keyheight = 0.3))

ggsave("./output/graphics/library-table.pdf", 
       width = 15 * 0.33, height = 7 * 0.33,
       device = cairo_pdf)
  























