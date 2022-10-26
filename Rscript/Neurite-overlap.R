
# Overlap with Neurite data -----------------------------------------------

library(tidyverse)
library(readxl)
library(patchwork)


# Prepare data ------------------------------------------------------------

## Mouse space ##

### Neurite data
neurite_enriched <- read_excel("./data/RNAseq/external_data/Supplementary online tables_neurite-enriched.xlsx",
                               skip = 1) %>% janitor::clean_names()

neurite_present <- read_excel("./data/RNAseq/external_data/Supplementary online tables_Chekulaeva_most-abundant.xlsx",
                              skip = 1) %>% janitor::clean_names()

### Glial data 
df <- readRDS("./data/summary_table.RDS")

### Summarise 
neurite_rna_present_n8 <- neurite_present %>% 
  filter(datasets_with_neurite_tpm_10 >= 8) %>% pull(gene_id)
neurite_translation_present_n2 <- neurite_present %>%
  filter(studies_with_neurite_ribosome_association >= 2) %>% pull(gene_id)
neurite_rna_enriched_n6 <- neurite_enriched %>%
  filter(datasets_with_significant_neurite_enrichment_p_0_1 >= 6) %>% pull(gene_id)

df_with_neurite <- df %>%
  mutate(neurite_rna_present = if_else(gene_id %in% neurite_rna_present_n8, TRUE, FALSE)) %>%
  mutate(neurite_translation_present = if_else(gene_id %in% neurite_translation_present_n2, TRUE, FALSE)) %>%
  mutate(neurite_rna_enriched = if_else(gene_id %in% neurite_rna_enriched_n6, TRUE, FALSE))

df_with_neurite_summary <- df_with_neurite %>%
  mutate(RNA_present = case_when(
    RNA_in_protrusion >= 8 & neurite_rna_present == TRUE ~ "Both",
    RNA_in_protrusion >= 8 ~ "Glia only",
    neurite_rna_present == TRUE ~ "Neurite only",
    TRUE ~ as.character(NA)
  )) %>%
  mutate(translation_present = case_when(
    translation_in_protrusion >= 4 & neurite_translation_present == TRUE ~ "Both",
    translation_in_protrusion >= 4 ~ "Glia only",
    neurite_translation_present == TRUE ~ "Neurite only",
    TRUE ~ as.character(NA)
  )) %>%
  mutate(RNA_enriched = case_when(
    enriched_in_protrusion >= 3 & neurite_rna_enriched == TRUE ~ "Both",
    enriched_in_protrusion >= 3 ~ "Glia only",
    neurite_rna_enriched == TRUE ~ "Neurite only",
    TRUE ~ as.character(NA)
  )) %>%
  dplyr::select(RNA_present, translation_present, RNA_enriched) %>%
  pivot_longer(cols = everything(),
               names_to = "type",
               values_to = "celltype") %>%
  group_by(type, celltype) %>%
  summarise(count = n()) %>% ungroup() %>%
  filter(!is.na(celltype)) %>%
  group_by(type) %>%
  mutate(total_count = sum(count)) %>% ungroup() %>%
  mutate(percentage = count/total_count * 100) %>%
  mutate(label = paste0(round(percentage, digits = 1), "%\n", "(", count, ")")) %>%
  mutate(celltype = str_replace_all(celltype, " ", "\n")) %>%
  mutate(type = str_replace_all(type, "_", "\n")) %>%
  mutate(type = str_replace_all(type, "trans", "Trans")) %>%
  mutate(type = fct_relevel(type, c("RNA\npresent", "RNA\nenriched", "Translation\npresent"))) %>%
  mutate(species = "Mouse")

## Drosophila space ##

diopt_conversion_table <- read_csv("./data/DIOPT/diopt_conversion_table.csv")

dmel_space_with_neurite <- diopt_conversion_table %>%
  left_join(df, by = c("mmus_gene_id" = "gene_id")) %>%
  mutate(neurite_rna_present = if_else(mmus_gene_id %in% neurite_rna_present_n8, TRUE, FALSE)) %>%
  mutate(neurite_translation_present = if_else(mmus_gene_id %in% neurite_translation_present_n2, TRUE, FALSE)) %>%
  mutate(neurite_rna_enriched = if_else(mmus_gene_id %in% neurite_rna_enriched_n6, TRUE, FALSE)) %>%
  group_by(dmel_gene_id, dmel_gene_name) %>%
  summarise(
    mmus_gene_ids = paste(mmus_gene_id, collapse = "|"),
    mmus_gene_names = paste(mmus_gene_name, collapse = "|"),
    RNA_in_protrusion = max(RNA_in_protrusion),
    enriched_in_protrusion = max(enriched_in_protrusion),
    translation_in_protrusion = max(translation_in_protrusion),
    enhanced_translation_in_protrusion = max(enhanced_translation_in_protrusion),
    neurite_rna_present = sum(neurite_rna_present),
    neurite_translation_present = sum(neurite_translation_present),
    neurite_rna_enriched = sum(neurite_rna_enriched)
  ) %>%
  ungroup()

# sanity check
# dmel_space %>%
#   filter(!is.na(RNA_in_protrusion)) %>%
#   pull(dmel_gene_name) %>% unique() %>% length()

dmel_space_with_neurite_summary <- dmel_space_with_neurite %>% 
  mutate(RNA_present = case_when(
    RNA_in_protrusion >= 8 & neurite_rna_present > 0 ~ "Both",
    RNA_in_protrusion >= 8 ~ "Glia only",
    neurite_rna_present > 0 ~ "Neurite only",
    TRUE ~ as.character(NA)
  )) %>%
  mutate(translation_present = case_when(
    translation_in_protrusion >= 4 & neurite_translation_present > 0 ~ "Both",
    translation_in_protrusion >= 4 ~ "Glia only",
    neurite_translation_present > 0 ~ "Neurite only",
    TRUE ~ as.character(NA)
  )) %>%
  mutate(RNA_enriched = case_when(
    enriched_in_protrusion >= 3 & neurite_rna_enriched > 0 ~ "Both",
    enriched_in_protrusion >= 3 ~ "Glia only",
    neurite_rna_enriched > 0 ~ "Neurite only",
    TRUE ~ as.character(NA)
  )) %>%
  dplyr::select(RNA_present, translation_present, RNA_enriched) %>%
  pivot_longer(cols = everything(),
               names_to = "type",
               values_to = "celltype") %>%
  group_by(type, celltype) %>%
  summarise(count = n()) %>% ungroup() %>%
  filter(!is.na(celltype)) %>%
  group_by(type) %>%
  mutate(total_count = sum(count)) %>% ungroup() %>%
  mutate(percentage = count/total_count * 100) %>%
  mutate(label = paste0(round(percentage, digits = 1), "%\n", "(", count, ")")) %>%
  mutate(celltype = str_replace_all(celltype, " ", "\n")) %>%
  mutate(type = str_replace_all(type, "_", "\n")) %>%
  mutate(type = str_replace_all(type, "trans", "Trans")) %>%
  mutate(type = fct_relevel(type, c("RNA\npresent", "RNA\nenriched", "Translation\npresent"))) %>%
  mutate(species = "Fly")

## Combine species ##

neurite_overlap_summary <- bind_rows(df_with_neurite_summary, dmel_space_with_neurite_summary)
View(neurite_overlap_summary)

# Plot --------------------------------------------------------------------
c("#E58606", "#5D69B1", "#52BCA3", "#99C945", "#CC61B0", "#24796C", "#DAA51B", "#2F8AC4", "#764E9F", "#ED645A", "#CC3A8E", "#A5AA99")

## RNA present only
fill_colour <- c("#E58606", "#52BCA3", "#2F8AC4")

neurite_overlap_summary %>%
  mutate(species = fct_rev(species)) %>%
  filter(type == "RNA\npresent") %>%
  ggplot(aes(
    x = species,
    y = count,
    fill = celltype
  )) +
  geom_col(width = 0.4, alpha = 1) +
  geom_text(aes(
    x = 1.4,
    label = celltype
  ), position = position_stack(vjust = 0.5), cex = 1.5) +
  geom_text(aes(
    x = 1.2,
    label = "-"
  ), position = position_stack(vjust = 0.5), cex = 1.5) +
  geom_label(aes(
    x = species,
    # y = count,
    label = label,
    group = celltype
  ), 
  label.padding = unit(0.1, "lines"),
  position = position_stack(vjust = 0.5), cex = 1.1, fill = "white", alpha = 0.5, fontface = "bold") + 
  labs(x = "",
       y = "Proportion of gene counts",
       fill = "") +
  scale_fill_manual(values = fill_colour) + 
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~ species, scales = "free") +
  theme_minimal(base_size = 6) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 6),
        panel.background = element_blank(),
        plot.background = element_blank())

ggsave("./output/graphics/overlap_with_neurite_data_bar.pdf",
       width = 5 * 0.33, height = 5.5 * 0.33, device = cairo_pdf)
  


## RNA present and enriched + translated
c("#E58606", "#5D69B1", "#52BCA3", "#99C945", "#CC61B0", "#24796C", "#DAA51B", "#2F8AC4", "#764E9F", "#ED645A", "#CC3A8E", "#A5AA99")

fill_colour <- c("#52BCA3", "#ED645A", "#2F8AC4")

neurite_overlap_summary %>%
  filter(species == "Mouse") %>%
  mutate(celltype = str_replace(celltype, "\n", " ")) %>%
  mutate(celltype = fct_relevel(celltype, c("Neurite only", "Both", "Glia only"))) %>%
  ggplot(aes(
    x = type,
    y = count,
    fill = celltype
  )) +
  geom_col(width = 0.7, alpha = 1) +
  geom_text(aes(
    x = 1.53,
    label = celltype
  ), position = position_stack(vjust = 0.5), cex = 1) +
  # geom_text(aes(
  #   x = 1.35,
  #   label = "|"
  # ), position = position_stack(vjust = 0.5), cex = 1) +
  geom_label(aes(
    x = type,
    # y = count,
    label = label,
    group = celltype
  ), 
  label.padding = unit(0.1, "lines"),
  position = position_stack(vjust = 0.5), label.size = 0.1, cex = 1.1, fill = "white", alpha = 0.5, 
  # fontface = "bold"
  ) + 
  labs(x = "",
       y = "Proportion of genes (Mouse)",
       fill = "") +
  scale_fill_manual(values = fill_colour) + 
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~ type, scales = "free", ncol = 1, strip.position = "left") +
  coord_flip() +
  theme_minimal(base_size = 6) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 4),
        panel.background = element_blank(),
        plot.background = element_blank()) -> mouse_plot

neurite_overlap_summary %>%
  filter(species == "Fly") %>%
  mutate(celltype = str_replace(celltype, "\n", " ")) %>%
    mutate(celltype = fct_relevel(celltype, c("Neurite only", "Both", "Glia only"))) %>%
  ggplot(aes(
    x = type,
    y = count,
    fill = celltype
  )) +
  geom_col(width = 0.7, alpha = 1) +
  geom_text(aes(
    x = 1.53,
    label = celltype
  ), position = position_stack(vjust = 0.5), cex = 1) +
  # geom_text(aes(
  #   x = 1.35,
  #   label = "|"
  # ), position = position_stack(vjust = 0.5), cex = 1) +
  geom_label(aes(
    x = type,
    # y = count,
    label = label,
    group = celltype
  ), 
  label.padding = unit(0.1, "lines"),
  position = position_stack(vjust = 0.5), label.size = 0.1, cex = 1.1, fill = "white", alpha = 0.5, 
  # fontface = "bold"
  ) + 
  labs(x = "",
       y = "Proportion of genes (Fly)",
       fill = "") +
  scale_fill_manual(values = fill_colour) + 
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~ type, scales = "free", ncol = 1, strip.position = "left") +
  coord_flip() +
  theme_minimal(base_size = 6) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 4),
        panel.background = element_blank(),
        plot.background = element_blank()) -> fly_plot

mouse_plot + fly_plot

ggsave("./output/graphics/overlap_with_neurite_data_full.pdf",
       width = 13 * 0.33, height = 4 * 0.33, device = cairo_pdf)



neurite_overlap_summary %>%
  

















