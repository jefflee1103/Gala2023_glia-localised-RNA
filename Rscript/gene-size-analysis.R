library(tidyverse)
library(valr)

gtf <- read_gtf("~/Documents/LCLIP/GENOMEDIR/Drosophila_melanogaster.BDGP6.28.99.chr.gtf")

genes <- filter(gtf, type == "gene" & gene_biotype == "protein_coding")




View(genes)










