library(tidyverse)
library(qs)

source("./Rscript/run_hypergeometric_test.R")

## Overlap between meta-analysis predicted and screen data
overlap <- 15

## Set size of predicted
set1 <- 1700

## Set size of screen data (dmel <-> mmus convertible)
gl_nd <- qread("./output/glia-localised-neurodegenerative-disease-incorporated.qs")
screen_paper <- read_csv("./data/screen-paper/screen_cpti-figureid-mapper.csv")

set2 <- filter(gl_nd, dmel_gene_id %in% screen_paper$FBgn_ID) %>% nrow()
set2

## Total background size (dmel <-> mmus convertible)
background <- nrow(gl_nd)
background

## Test enrichment
run_hypergeometric_test(
    11, # Overlap between prediction & screen glial localised genes
    1700, # meta-analyisis predicted genes
    15, # screen glial localised genes (mouse<->fly convertible)
    7011 # mouse<->fly convertible homologs
)

