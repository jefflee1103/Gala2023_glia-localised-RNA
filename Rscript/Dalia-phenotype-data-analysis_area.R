##
## Analyse Dalia's phenotype data
##

library(tidyverse)
library(rstatix)
library(patchwork)
library(furrr)
library(colorspace)
plan(multisession, workers = 4)

## Genes to test 
goi <- c("lac", "Pdi", "Gs2", "nrv2", "alpha-Cat", "ATP-alpha", "nrg", "shot", "Vha55", "flo2")

# ----- Developmental defects in glial/neuronal area 

## 
dev_area_raw <- list.files(
    "./data/Dalias_Data/Functional_data_prediction_paper/Developmental_areas_analysis",
    pattern = "*.csv",
    full.names = TRUE
  ) %>%
  map(read_csv)

## 
dev_area_statistics <- dev_area_raw %>%
    map_dfr(~ {
        dev_area_cleaned <- .x %>%
            mutate(Label = str_replace(Label, "_C[:digit:].tif", "")) %>%
            pivot_wider(id_cols = c(Label, condition), values_from = sum_Area, names_from = channel) %>%
            dplyr::rename("Neurite" = HRP) %>%
            mutate(ratio_Glia_Neurite = Glia / Neurite) %>%
            pivot_longer(cols = c(Glia, Neurite, ratio_Glia_Neurite), names_to = "type")

        dev_area_mean <- dev_area_cleaned %>%
            group_by(condition, type) %>%
            summarise(
                mean = mean(value)
            ) %>%
            ungroup()

        dev_area_foldchange <- left_join(
            filter(dev_area_mean, condition != "Control"),
            filter(dev_area_mean, condition == "Control") %>%
                dplyr::select(type, "control_mean" = mean),
            by = "type"
        ) %>%
            mutate(foldchange = mean / control_mean)

        dev_area_statistics <- dev_area_cleaned %>%
            group_by(type) %>%
            t_test(value ~ condition, ref.group = "Control") %>%
            dplyr::select(type, "condition" = group2, statistic, p) %>%
            left_join(dev_area_foldchange, by = c("type", "condition"))
    })

##
dev_area_statistics_for_plotting <- dev_area_statistics %>%
  mutate(p = case_when(
    p < 0.0001 ~ "**** p < 0.0001",
    p < 0.001 ~ "*** p < 0.001",
    p < 0.01 ~ "** p < 0.01",
    p < 0.05 ~ "* p < 0.05",
    TRUE ~ "N.S."
  )) %>%
  mutate(type = case_when(
    type == "Glia" ~ "Glial area",
    type == "Neurite" ~ "Neurite area",
    type == "ratio_Glia_Neurite" ~ "Ratio Glial/Neurite"
  )) %>%
  mutate(condition = str_replace(condition, "_VDRC", ""))

colours <- c(sequential_hcl(n = 6, palette = "Inferno", rev = FALSE)[3:5], "gray70")
gene_order <- dev_area_statistics_for_plotting %>%
  filter(type == "Glial area") %>%
  arrange(desc(foldchange)) %>% pull(condition)

dev_area_statistics_for_plotting %>%
  mutate(condition = fct_relevel(condition, gene_order)) %>%
  ggplot(aes(x = condition, y = foldchange, colour = p)) +
  geom_segment(aes(xend = condition, y = 1, yend = foldchange), colour = "gray50", size = 0.3) + 
  geom_point(size = 1.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", colour = "gray70", size = 0.2) +
  coord_flip(ylim = c(0.5, 1.7)) +
  labs(
    x = "", 
    y = "Average foldchange over Control",
    colour = "P-value"
  ) + 
  scale_colour_manual(values = colours) + 
  facet_wrap(~ type) +
  theme_bw(base_size = 6) + 
  theme(
    legend.text = element_text(size = 4),
    legend.title = element_text(size = 5),
    axis.text.x = element_text(size = 3)
    ) +
  guides(colour = guide_legend(keywidth = 0.1, keyheight = 0.5))

ggsave("./output/graphics/phenotype_area_plot_dev.pdf", 
       width = 11 * 0.3937, height = 5 * 0.38937,
       device = cairo_pdf)


# ----- Post-kstim defects in glial/neuronal area 

## 
kstim_area_raw <- list.files(
    "./data/Dalias_Data/Functional_data_prediction_paper/Kstim_all_genes_areas_analysis/",
    pattern = "*.csv",
    full.names = TRUE
  ) %>%
  map(read_csv)

## 
kstim_area_statistics <- kstim_area_raw %>%
    map_dfr(~ {
        kstim_area_cleaned <- .x %>%
            mutate(Label = str_replace(Label, "_C[:digit:].tif", "")) %>%
            pivot_wider(id_cols = c(Label, condition), values_from = sum_Area, names_from = channel) %>%
            dplyr::rename("Neurite" = HRP) %>%
            mutate(ratio_Glia_Neurite = Glia / Neurite) %>%
            pivot_longer(cols = c(Glia, Neurite, ratio_Glia_Neurite), names_to = "type")

        kstim_area_mean <- kstim_area_cleaned %>%
            group_by(condition, type) %>%
            summarise(
                mean = mean(value)
            ) %>%
            ungroup()

        kstim_area_foldchange <- left_join(
            filter(kstim_area_mean, condition != "Control"),
            filter(kstim_area_mean, condition == "Control") %>%
                dplyr::select(type, "control_mean" = mean),
            by = "type"
        ) %>%
            mutate(foldchange = mean / control_mean)

        kstim_area_statistics <- kstim_area_cleaned %>%
            group_by(type) %>%
            t_test(value ~ condition, ref.group = "Control") %>%
            dplyr::select(type, "condition" = group2, statistic, p) %>%
            left_join(kstim_area_foldchange, by = c("type", "condition"))
    })

## Replace if multi-replicate data is available
ksim_area_raw_triplicate <- list.files(
    "./data/Dalias_Data/Functional_data_prediction_paper/KStim_triplicates_areas_analysis",
    pattern = "*.csv",
    full.names = TRUE
) %>%
  map(read_csv)

kstim_area_statistics_triplicate <- ksim_area_raw_triplicate %>%
    map_dfr(~ {
        kstim_area_cleaned <- .x %>%
            mutate(Label = str_replace(Label, "_C[:digit:].tif", "")) %>%
            pivot_wider(id_cols = c(Label, condition), values_from = sum_Area, names_from = channel) %>%
            dplyr::rename("Neurite" = HRP) %>%
            mutate(ratio_Glia_Neurite = Glia / Neurite) %>%
            pivot_longer(cols = c(Glia, Neurite, ratio_Glia_Neurite), names_to = "type")

        kstim_area_mean <- kstim_area_cleaned %>%
            group_by(condition, type) %>%
            summarise(
                mean = mean(value)
            ) %>%
            ungroup()

        kstim_area_foldchange <- left_join(
            filter(kstim_area_mean, condition != "Control"),
            filter(kstim_area_mean, condition == "Control") %>%
                dplyr::select(type, "control_mean" = mean),
            by = "type"
        ) %>%
            mutate(foldchange = mean / control_mean)

        kstim_area_statistics <- kstim_area_cleaned %>%
            group_by(type) %>%
            t_test(value ~ condition, ref.group = "Control") %>%
            dplyr::select(type, "condition" = group2, statistic, p) %>%
            left_join(kstim_area_foldchange, by = c("type", "condition"))
    })

### Remove CG1648 and combine with the summary data 
kstim_area_statistics_updated <- kstim_area_statistics_triplicate %>%
  filter(!str_detect(condition, "CG")) %>%
  bind_rows(
    filter(kstim_area_statistics, !(condition %in% c("lac", "gs2", "Pdi")))
  )

## Create a df for plotting - clean factors 
kstim_area_statistics_for_plotting <- kstim_area_statistics_updated %>%
  mutate(p = case_when(
    p < 0.0001 ~ "**** p < 0.0001",
    p < 0.001 ~ "*** p < 0.001",
    p < 0.01 ~ "** p < 0.01",
    p < 0.05 ~ "* p < 0.05",
    TRUE ~ "N.S."
  )) %>%
  mutate(type = case_when(
    type == "Glia" ~ "Glial area",
    type == "Neurite" ~ "Neurite area",
    type == "ratio_Glia_Neurite" ~ "Ratio Glial/Neurite"
  ))

colours <- c(sequential_hcl(n = 6, palette = "Inferno", rev = FALSE)[3:5], "gray70")
gene_order <- kstim_area_statistics_for_plotting %>%
  filter(type == "Glial area") %>%
  arrange(desc(foldchange)) %>% pull(condition)

kstim_area_statistics_for_plotting %>%
  mutate(condition = fct_relevel(condition, gene_order)) %>%
  ggplot(aes(x = condition, y = foldchange, colour = p)) +
  geom_segment(aes(xend = condition, y = 1, yend = foldchange), colour = "gray50", size = 0.3) + 
  geom_point(size = 1.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", colour = "gray70", size = 0.2) +
  coord_flip(
    # ylim = c(0.5, 1.7)
    ) +
  labs(
    x = "", 
    y = "Average foldchange over Control",
    colour = "P-value"
  ) + 
  scale_colour_manual(values = colours) + 
  facet_wrap(~ type) +
  theme_bw(base_size = 6) + 
  theme(
    legend.text = element_text(size = 4),
    legend.title = element_text(size = 5),
    axis.text.x = element_text(size = 3)
    ) +
  guides(colour = guide_legend(keywidth = 0.1, keyheight = 0.5))

ggsave("./output/graphics/phenotype_area_plot_kstim.pdf", 
       width = 11 * 0.3937, height = 5 * 0.38937,
       device = cairo_pdf)



