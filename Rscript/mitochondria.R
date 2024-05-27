# Intersect 1700 genes with mitochondrial function

# Environment ------ 
library(tidyverse)
library(rstatix)
library(qs)
library(colorspace)

# GO terms ----- 
mito_go <- list.files("./data/mito/", pattern = "*.txt", full.names = TRUE) %>%
  set_names(~ basename(tools::file_path_sans_ext(.x))) %>%
  map_dfr(~{
    read_tsv(.x, col_names = FALSE)
  }, .id = "filename") %>%
  dplyr::rename("gene_id" = X1) %>%
  separate_wider_delim(cols = "filename", delim = "_", names = c("go_id", "go_name"), cols_remove = TRUE)

# Annotate with 1700 -----
gl_nd <- qread("./output/glia-localised-neurodegenerative-disease-incorporated.qs")
id_interest <- gl_nd %>%
  filter(rna_in_protrusion >= 8) %>%
  filter(`subperineurial glial cell` == TRUE | 
           `perineurial glial sheath` == TRUE |
           `adult brain perineurial glial cell` == TRUE |
           `ensheathing glial cell` == TRUE)

mito_go_anno <- mito_go %>%
  left_join(id_interest, by = c("gene_id" = "dmel_gene_id")) %>%
  mutate(is_1700 = if_else(is.na(dmel_gene_name), FALSE, TRUE))

# Prepare csv of list of genes -----
output_csv <- mito_go_anno %>%
  filter(!is.na(dmel_gene_name)) %>%
  group_by(go_name) %>%
  summarise(
    dmel_gene_name = str_c(sort(unique(dmel_gene_name)), collapse = ", "),
    mmus_gene_name = str_c(sort(unique(mmus_gene_symbols)), collapse = ", "),
    )

write_csv(output_csv, "./output/mito_function_transcripts_localised_to_glialprotrusion.csv")

# Summary with convertible genes -----
## phyper(overlap - 1, group2, Total - group2, group1, lower.tail = FALSE)
mito_go_anno_summary <- mito_go_anno %>%
  filter(gene_id %in% gl_nd$dmel_gene_id) %>%
  group_by(go_name) %>%
  summarise(
    go_gene_total_count = n(),
    go_in1700_count = sum(is_1700)
  ) %>%
  mutate(hyper_p = phyper(
    go_in1700_count - 1,
    go_gene_total_count,
    nrow(gl_nd) - go_gene_total_count,
    nrow(id_interest),
    lower.tail = FALSE
  )) %>%
  adjust_pvalue(p.col = "hyper_p", output.col = "hyper_padj", method = "bonferroni") %>%
  add_significance(p.col = "hyper_padj") %>%
  mutate(foldchange =
    (go_in1700_count / nrow(id_interest)) / 
      (go_gene_total_count / nrow(gl_nd))
  ) %>%
  mutate(proportion = str_c(as.character(go_in1700_count), "/", as.character(go_gene_total_count))) %>%
  arrange(hyper_padj)

go_order <- mito_go_anno_summary$go_name
colours <- sequential_hcl(palette = "inferno", n = 6)[c(3,4,5)]

mito_go_anno_summary %>%
  mutate(go_name = fct_relevel(go_name, rev(go_order))) %>%
  ggplot(aes(x = foldchange, y = go_name, fill = hyper_padj.signif)) +
  geom_col(width = 0.5) + 
  geom_text(aes(label = proportion, x = foldchange - 0.2), cex = 3) + 
  labs(title = "Enriched mitochondrial GO terms", 
       subtitle = "1700 glial protrusion genes",
       x = "Fold enrichment", y = "", fill = "Overlap padj") +
  scale_fill_manual(values = c(colours, "gray80")) + 
  theme_minimal()

ggsave("./output/graphics/mito_go_barplot.pdf", width = 7, height = 5)

# Network plot -----
library(ggraph)
library(tidygraph)
library(ggrepel)

library(enrichplot)

clusterproflier_singles <- qread("~/Documents/LCLIP/analysis_data/go/clusterprofiler_singles.qs")

x <- clusterproflier_singles$bp$L1_Imp

inputList <- x

ldf <- lapply(seq_len(length(inputList)), function(i) {
  data.frame(categoryID=rep(names(inputList[i]),
                            length(inputList[[i]])),
             Gene=inputList[[i]])
})

do.call('rbind', ldf)

extract_geneSets <- function(x, n) {
  n <- update_n(x, n)
  geneSets <- geneInCategory(x) ## use core gene for gsea result
  y <- as.data.frame(x)
  geneSets <- geneSets[y$ID]
  names(geneSets) <- y$Description
  if (is.numeric(n)) {
    return(geneSets[1:n])
  }
  return(geneSets[n]) ## if n is a vector of Description
}

update_n <- function(x, showCategory) {
  if (!is.numeric(showCategory)) {
    return(showCategory)
  }
  
  ## geneSets <- geneInCategory(x) ## use core gene for gsea result
  n <- showCategory
  if (nrow(x) < n) {
    n <- nrow(x)
  }
  
  return(n)
}
list2graph <- function(inputList) {
  x <- list2df(inputList)
  g <- graph.data.frame(x, directed=FALSE)
  return(g)
}


list2df <- function(inputList) {
  # ldf <- lapply(1:length(inputList), function(i) {
  ldf <- lapply(seq_len(length(inputList)), function(i) {
    data.frame(categoryID=rep(names(inputList[i]),
                              length(inputList[[i]])),
               Gene=inputList[[i]])
  })
  
  do.call('rbind', ldf)
}
highschool %>% head()
graph <- as_tbl_graph(highschool) 

graph[[2]]

a <- list2df(geneSets) %>% 
  as_tibble() %>%
  as_tbl_graph()

ggraph(a, layout = "kk") +
  geom_edge_link() + 
  geom_node_point(aes_(color=~I("#E5C494")))





showCategory <- 5

extract_geneSets(x, n = 5)

geneInCategory(x)
as.data.frame(x) -> df
View(df)



names(geneSets)

geneSets <- setNames(strsplit(as.character(y_union$geneID), "/",
                              fixed = TRUE), y_union$Description)

create_notable('bull') |>
  mutate(class = sample(letters[1:3], n(), replace = TRUE))

## 
mito_graph <- mito_go_anno %>%
  mutate(go_name = if_else(str_detect(go_name, "mitochondrial ATP synthesis"), 
                           "mitochondrial ATP synthesis", 
                           go_name)) %>%
  filter(!is.na(dmel_gene_name)) %>%
  dplyr::select(go_name, dmel_gene_name) %>%
  as_tbl_graph(directed = FALSE)

go_categories <- unique(mito_go_anno$go_name)

p <- mito_graph %>%
  ggraph(layout = "kk") +
  geom_edge_link(alpha = 0.1) +
  theme_void()

p + 
  geom_node_point(
    data = p$data[-c(1:length(go_categories)),],
    colour = "#d18975",
    cex = 2.5
  ) + 
  geom_node_label(
    data = p$data[1:length(go_categories),],
    aes(label = name),
    cex = 3.5
  ) +
  geom_node_text(
    data = p$data[-c(1:length(go_categories)),],
    aes(label = name),
    repel = TRUE, cex = 3
  ) + 
  labs(title = "Glial localised transcripts involved in mitochondrial organisation/function") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave("./output/graphics/mito_go_networkplot.pdf", width = 18, height = 10)



