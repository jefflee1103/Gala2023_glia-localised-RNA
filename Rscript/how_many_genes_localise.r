library(tidyverse)
library(qs)

df <- readRDS("./data/RNAseq/quant_results/all-glia_avgTPM_wide.RDS")

df_long <- df %>%
  pivot_longer(cols = !gene_id, names_to = "library", values_to = "TPM")

df_long_tidy <- df_long %>%
  mutate(glia_type = str_replace(library, "_protrusion", "")) %>%
  mutate(glia_type = str_replace(glia_type, "_soma", "")) %>%
  mutate(compartment = if_else(str_detect(library, "protrusion"), "protrusion", "soma"))

count_table <- df_long_tidy$glia_type %>%
  unique() %>%
  set_names() %>%
  map_dfr(~{
    all_sequenceable <- df_long_tidy %>%
      filter(glia_type == .x) %>%
      filter(TPM > 0) %>%
      pull(gene_id) %>% unique() %>% length()

    expressed_gene_count <- df_long_tidy %>%
      filter(glia_type == .x) %>%
      filter(TPM >= 1) %>%
      pull(gene_id) %>% unique() %>% length()
    
    protrusion_expressed_gene_count <- df_long_tidy %>%
      filter(glia_type == .x) %>%
      filter(compartment == "protrusion") %>%
      filter(TPM >= 1) %>%
      pull(gene_id) %>% unique() %>% length()
    
    compartment_number <- df_long_tidy %>%
      filter(glia_type == .x) %>%
      pull(compartment) %>% unique() %>% length()

    tibble(
        sequenceable = all_sequenceable,
        expressed = expressed_gene_count,
        protrusion = protrusion_expressed_gene_count,
        compartments = compartment_number
    )
  }, .id = "glia_type", .progress = TRUE)

pct_table <- count_table %>%
  mutate(percentage_protrusion = case_when(
    compartments == 2 ~ (protrusion / expressed) * 100,
    compartments == 1 ~ (protrusion / sequenceable) * 100
  ))



get_pct_table <- function(cutoff){
  count_table <- df_long_tidy$glia_type %>%
    unique() %>%
    set_names() %>%
    map_dfr(~{
      all_sequenceable <- df_long_tidy %>%
        filter(glia_type == .x) %>%
        filter(TPM > 0) %>%
        pull(gene_id) %>% unique() %>% length()
      
      expressed_gene_count <- df_long_tidy %>%
        filter(glia_type == .x) %>%
        filter(TPM >= cutoff) %>%
        pull(gene_id) %>% unique() %>% length()
      
      protrusion_expressed_gene_count <- df_long_tidy %>%
        filter(glia_type == .x) %>%
        filter(compartment == "protrusion") %>%
        filter(TPM >= cutoff) %>%
        pull(gene_id) %>% unique() %>% length()
      
      compartment_number <- df_long_tidy %>%
        filter(glia_type == .x) %>%
        pull(compartment) %>% unique() %>% length()
      
      tibble(
        sequenceable = all_sequenceable,
        expressed = expressed_gene_count,
        protrusion = protrusion_expressed_gene_count,
        compartments = compartment_number
      )
    }, .id = "glia_type", .progress = TRUE)
  
  pct_table <- count_table %>%
    mutate(percentage_protrusion = case_when(
      compartments == 2 ~ (protrusion / expressed) * 100,
      compartments == 1 ~ (protrusion / sequenceable) * 100
    ))
  
  return(pct_table)
}

get_pct_table_range <- function(cutoff_start, cutoff_end){
  count_table <- df_long_tidy$glia_type %>%
    unique() %>%
    set_names() %>%
    map_dfr(~{
      all_sequenceable <- df_long_tidy %>%
        filter(glia_type == .x) %>%
        filter(TPM > 0) %>%
        pull(gene_id) %>% unique() %>% length()
      
      expressed_gene_count <- df_long_tidy %>%
        filter(glia_type == .x) %>%
        filter(TPM >= cutoff_start, TPM <= cutoff_end) %>%
        pull(gene_id) %>% unique() %>% length()
      
      protrusion_expressed_gene_count <- df_long_tidy %>%
        filter(glia_type == .x) %>%
        filter(compartment == "protrusion") %>%
        filter(TPM >= cutoff_start, TPM <= cutoff_end) %>%
        pull(gene_id) %>% unique() %>% length()
      
      compartment_number <- df_long_tidy %>%
        filter(glia_type == .x) %>%
        pull(compartment) %>% unique() %>% length()
      
      tibble(
        sequenceable = all_sequenceable,
        expressed = expressed_gene_count,
        protrusion = protrusion_expressed_gene_count,
        compartments = compartment_number
      )
    }, .id = "glia_type", .progress = TRUE)
  
  pct_table <- count_table %>%
    mutate(percentage_protrusion = case_when(
      compartments == 2 ~ (protrusion / expressed) * 100,
      compartments == 1 ~ (protrusion / sequenceable) * 100
    ))
  
  return(pct_table)
}

quantiles_vector <- df_long_tidy %>%
  filter(TPM > 0) %>%
  pull(TPM) %>%
  log10() %>% 
  quantile() %>%
  as.vector()

quantile_range_list <- list()
for (i in 1:(length(quantiles_vector) - 1)){
  start <- 10^(quantiles_vector[i])
  end <- 10^(quantiles_vector[i+1])
  output <- get_pct_table_range(cutoff_start = start, cutoff_end = end) %>%
    mutate(quantile = i)
  quantile_range_list[[i]] <- output
}
quantile_range_df <- purrr::reduce(quantile_range_list, bind_rows)

quantile_range_df %>%
  ggplot(aes(x = as.factor(quantile), y = percentage_protrusion, colour = as.factor(quantile))) +
  geom_boxplot(width = 0.5) + 
  geom_jitter(colour = "black", position = position_jitter(width = 0.2)) +
  labs(title = "No correlation between expression level\nand protrusion localisation",
       subtitle = "",
       x = "TPM quantiles",
       y = "% localised transcripts"
       ) + 
  coord_cartesian(ylim = c(0, 100)) + 
  theme_classic(base_size = 17) +
  theme(legend.position = "none")





pct_table_list <- df_long_tidy %>%
  filter(TPM > 0) %>%
  pull(TPM) %>%
  log10() %>% 
  quantile() %>%
  as.list() %>%
  map_dfr(~{
    get_pct_table(cutoff = 10^(.x))
  }, .id = "quantiles")
  

pct_table_list %>%
  group_by(quantiles) %>%
  summarise(mean_localisation_pct = mean(percentage_protrusion)) %>%
  arrange(quantiles)




quantile(log10(df_long_tidy$TPM + 1))







df_long$library %>% unique()



View(df_long)







