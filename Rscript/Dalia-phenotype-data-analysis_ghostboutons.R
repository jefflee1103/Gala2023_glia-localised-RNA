##
## Analyse Dalia's phenotype data
##

library(tidyverse)
library(rstatix)
library(patchwork)
library(furrr)
library(colorspace)
plan(multisession, workers = 4)

# ----- Kstim ghost bouton counts 

## Helper function
get_summary_results <- function(list_of_df) {
  map(list_of_df, ~ {
    kstim_mean <- .x %>%
      group_by(condition) %>%
      summarise(mean_bouton_count = mean(bouton_count))

    kstim_foldchange <- bind_cols(
      filter(kstim_mean, condition != "Control"),
      filter(kstim_mean, condition == "Control") %>%
        dplyr::select("control_mean" = mean_bouton_count)
    ) %>%
      mutate(foldchange = mean_bouton_count / control_mean)

    kstim_statistics_wilcox <- .x %>%
      wilcox_test(bouton_count ~ condition, ref.group = "Control") %>%
      add_significance() %>%
      dplyr::select("condition" = group2, statistic, p, starts_with("n")) %>%
      left_join(kstim_foldchange, by = "condition")
    
    kstim_statistics_ttest <- .x %>%
      t_test(bouton_count ~ condition, ref.group = "Control") %>%
      add_significance() %>%
      dplyr::select("condition" = group2, statistic, p) %>%
      left_join(kstim_foldchange, by = "condition")


    output <- list(kstim_statistics_wilcox, kstim_statistics_ttest, .x) %>%
      set_names(c("Wilcox", "t-test", "full_df"))
  })
}

## Genes with single repilicates
kstim_raw <- list.files(
    "./data/Dalias_Data/Functional_data_prediction_paper/KStim_all_genes_ghost_boutons_analysis",
    pattern = "*.csv",
    full.names = TRUE
  ) %>%
  set_names() %>%
  map(read_csv) %>%
  imap(~ dplyr::select(.x, c(Set, Segment, contains("NMJ"), Boutons)) %>%
         setNames(c("condition", "larval_segment", "nmj_number", "bouton_count")) %>%
         mutate(condition = case_when(
          condition == "Control" ~ "Control",
          str_detect(condition, "RNAi") ~ str_replace(basename(.y), ".csv", ""),
          TRUE ~ condition
         ))
         )
  
kstim_statistics <- get_summary_results(kstim_raw) 

## Create summary dataframe
kstim_wilcox_df <- kstim_statistics %>%
  map_dfr(~ pluck(.x, "Wilcox")) %>%
  filter(!str_detect(condition, "CG")) %>%
  filter(condition != "cold")


## Correct gene names
t2g <- read_csv("./data/Flybase/Dmel_tx2gene_ENSEMBL_v99.csv")
id_change <- tibble(
  input = kstim_wilcox_df$condition,
  ens99_name = c("alpha-Cat", "Gs2", "Lac", "Pdi", "Atpalpha", "Flo2", "Nrg", "nrv2", "shot", "Vha55")
)

## 
df_for_plotting <- kstim_wilcox_df %>%
  left_join(id_change, by = c("condition" = "input")) %>%
  mutate(p_legend = case_when(
    p < 0.0001 ~ "**** < 0.0001",
    p < 0.001 ~ "*** < 0.001",
    p < 0.01 ~ "** < 0.01",
    p < 0.05 ~ "* < 0.05",
    TRUE ~ "N.S."
  )) %>%
  mutate(log2foldchange = log2(foldchange)) %>%
  mutate(p_label = str_c("p=", signif(p, digits = 2))) %>%
  mutate(p_label_pos = case_when(
    log2foldchange >= 0 ~ log2foldchange + 0.08,
    log2foldchange < 0 ~ log2foldchange - 0.08
  ))

##
colours <- c(sequential_hcl(n = 5, palette = "Inferno", rev = FALSE)[4], "gray70")

df_for_plotting %>%
  ggplot(aes(x = reorder(ens99_name, -foldchange), y = log2foldchange, fill = p_legend)) + 
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, colour = "gray20", size = 0.1) +
  geom_text(aes(y = p_label_pos, label = p_label), cex = 1.4) +
  coord_cartesian(ylim = c(-1.8, 1.7)) +
  labs(
    x = "", 
    y = "log2 Foldchange over Control",
    fill = "P-value"
  ) + 
  scale_fill_manual(values = colours) +
  theme_bw(base_size = 7) + 
  theme(
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    axis.text.x = element_text(size = 5)
    ) +
  guides(fill = guide_legend(keywidth = 0.3, keyheight = 0.3))

ggsave("./output/graphics/phenotype_boutons_plot_kstim.pdf", 
       width = 10.2 * 0.3937, height = 5.5 * 0.38937,
       device = cairo_pdf)

##
kstim_normalised_bouton_df <- kstim_full_df %>% 
  filter(condition != "Control") %>%
  left_join(
    dplyr::select(kstim_wilcox_df, c(condition, p , control_mean, mean_bouton_count, foldchange)),
    by = "condition"
  ) %>%
  mutate(normalised_bouton_count = (bouton_count + 1) / (control_mean)) %>%
  mutate(p_legend = case_when(
    p < 0.0001 ~ "**** < 0.0001",
    p < 0.001 ~ "*** < 0.001",
    p < 0.01 ~ "** < 0.01",
    p < 0.05 ~ "* < 0.05",
    TRUE ~ "N.S."
  )) %>%
  mutate(p_label = str_c("p=", signif(p, digits = 2))) 

kstim_normalised_bouton_df %>%
  ggplot(aes(x = reorder(condition, -foldchange), y = log2(normalised_bouton_count), fill = p_legend)) +
  geom_bar(stat = "summary", width = 0.6) +
  geom_linerange(stat = "summary", size = 0.25) +
  coord_cartesian(ylim = c(-2, 1.7)) +
  labs(
    x = "", 
    y = "log2 Bounton count foldchange over Control",
    fill = "P-value"
  ) + 
  scale_fill_manual(values = colours) +
  theme_bw(base_size = 6) + 
  theme(
    legend.text = element_text(size = 4),
    legend.title = element_text(size = 5)
    ) +
  guides(fill = guide_legend(keywidth = 0.3, keyheight = 0.3))

ggsave("./output/graphics/phenotype_boutons_plot_kstim_with_errorbar.pdf", 
       width = 10 * 0.3937, height = 5.5 * 0.38937,
       device = cairo_pdf)



















