library(SCopeLoomR)
library(tidyverse)
library(qs)
library(patchwork)
library(colorspace)

## Get glia FCA data and fiter its annotation embeddings
loom <- open_loom('./data/FCA/s_fca_biohub_glial_cell_10x.loom')
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
glia_anno <- c(
  "perineurial glial sheath",
  "subperineurial glial cell",
  "ensheathing glial cell",
  "adult antenna glial cell",
  "adult brain cell body glial cell",
  "adult brain perineurial glial cell",
  "adult glial cell",
  "adult lamina epithelial/marginal glial cell",
  "adult optic chiasma glial cell",
  "adult reticular neuropil associated glial cell",
  "cell body glial cell",
  "CNS surface associated glial cell",
  "optic-lobe-associated cortex glial cell",
  "peripheral glial cell"
)
all_glia_cell_anno <- filter(cell_anno, annotation %in% glia_anno)

## Plot (3 + 11 cell types)

colours <- c(sequential_hcl(3, palette = "Peach"), sequential_hcl(11, palette = "Teal"))
colours <- c(sequential_hcl(6, palette = "OrRd")[1:3], sequential_hcl(11, palette = "TealGrn"))

all_glia_cell_anno %>%
  mutate(annotation = fct_relevel(annotation, glia_anno)) %>%
  arrange(match(annotation, rev(glia_anno))) %>%
  ggplot(aes(x = hvg_t_sne_x,
             y = hvg_t_sne_y,
             colour = annotation)) +
  labs(colour = "",
    x = "t-SNE x",
    y = "t-SNE y") + 
  geom_point(alpha = 0.5, size = 0.001, stroke = 0) +
  scale_colour_manual(values = colours) +
  theme_classic(base_size = 6) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 3)
    ) +
  guides(colour = guide_legend(
    keywidth = 0.3,
    keyheight = 0.3,
    override.aes = list(size = 1, alpha = 0.8)
  ))

ggsave("./output/graphics/FCA-plot.pdf",
       width = 8 * 0.33, height = 5 * 0.33, device = cairo_pdf)
 




hcl_palettes(type = "quantitative")
