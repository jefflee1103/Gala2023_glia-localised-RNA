library(tidyverse)
library(patchwork)
library(qs)
library(colorspace)

## Mouse summary table from Dalia's RNA synaptic glia repo
rna_present <- readRDS("./data/summary_table.RDS") %>% 
  janitor::clean_names() %>%
  filter(rna_in_protrusion >= 8) 


## DIOPT data
diopt_conversion_table <- read_csv("./data/DIOPT/diopt_conversion_table.csv")
diopt_raw <- qread("./data/DIOPT/mouse96-to-fly/DIOPT_mouse-to-fly.qs")

## Get diopt summary
dmel_rna_present <- rna_present %>%
  filter(rna_in_protrusion >= 8) %>%
  left_join(diopt_raw, by = c("gene_id" = "mmus_gene_id"))

dmel_rna_present_diopt_summary <- dmel_rna_present %>%
  group_by(dmel_gene_id, dmel_gene_name) %>%
  summarise(max_diopt_score = max(diopt_score)) %>%
  ungroup() %>%
  group_by(max_diopt_score) %>%
  summarise(count = n()) %>%
  mutate(cumsum = cumsum(count)) %>%
  mutate(diopt_label = if_else(max_diopt_score %in% c(1, 8, 18), as.character(max_diopt_score), NA_character_)) %>%
  mutate(dash = if_else(!is.na(diopt_label), "-", NA_character_))

##
colours <- c(sequential_hcl(20, palette = "OrRd")[5:15], sequential_hcl(12, palette = "TealGrn")[1:7])

tibble(x = "x", y = nrow(rna_present)) %>%
  ggplot(aes(x = x, y = y)) +
  geom_col(width = 0.45, fill = "#E58606") + 
  coord_cartesian(ylim = c(0, 6000)) +
  labs(title = "Mouse", x = "", y = "Gene count") + 
  theme_minimal(base_size = 6) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 3),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank()) -> p1


dmel_rna_present_diopt_summary %>%
  filter(!is.na(max_diopt_score)) %>%
  mutate(max_diopt_score = as.factor(max_diopt_score)) %>%
  mutate(max_diopt_score = fct_rev(max_diopt_score)) %>%
  ggplot(aes(x = "", y = count, fill = as.factor(max_diopt_score))) +
  geom_col(position = "stack", width = 0.45) +
  geom_text(aes(
    x = 1.35,
    y = cumsum,
    label = diopt_label
  ), cex = 1.5) +
  geom_text(aes(
    x = 1.23,
    y = cumsum,
    label =  dash
  ), cex = 1.5) +
  scale_fill_manual(values = colours) +
  labs(title = "Fly", 
    x = "",
    y = "") +
  coord_cartesian(ylim = c(0, 6000)) +
  theme_minimal(base_size = 6) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank()) -> p2

p1 + p2

ggsave("./output/graphics/diopt-coversion-plot.pdf",
       width = 6 * 0.33, height = 5.5 * 0.33, device = cairo_pdf)

nrow(rna_present)




