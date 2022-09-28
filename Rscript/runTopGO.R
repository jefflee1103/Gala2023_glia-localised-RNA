# ------------------------------------------------------------------------
# Performing GO analysis with 'TopGo' package.
# Revised: 2020.09.01
# ------------------------------------------------------------------------

### Environment 
library(topGO)
library(tidyverse)
# need to load appropriate species database
library(org.Dm.eg.db)


### Set runTopGO fuction 
runTopGO <- function(ontology, species, backgroundGenes, targetedGenes) 
{
  category <- data.frame(geneId = backgroundGenes)
  category$targeted <- 0
  category[category$geneId %in% targetedGenes, ]$targeted <- 1
  category$geneId <- gsub("\\..*$", "", category$geneId)
  allGenes <- category$targeted
  names(allGenes) <- category$geneId
  if (species == "human") {
    speciesOrgdb <- "org.Hs.eg.db"
  }
  else if (species == "mouse") {
    speciesOrgdb <- "org.Mm.eg.db"
  }
  else if (species == "fly") {
    speciesOrgdb <- "org.Dm.eg.db"
  }
  else if (species == "worm") {
    speciesOrgdb <- "org.Ce.eg.db"
  }
  else {
    stop("Cannot do GO term analysis for species except:\n          human, worm, fly, and mouse\n")
  }
  GOdata <- new("topGOdata", ontology = ontology, allGenes = allGenes, 
                geneSelectionFun = function(x) {
                  return(x == 1)
                }, nodeSize = 20, description = "Test", annot = annFUN.org, 
                mapping = speciesOrgdb, ID = "Ensembl")
  resultFisher <- topGO::runTest(object = GOdata, algorithm = "classic", 
                                 statistic = "fisher")
  goResults <- topGO::GenTable(object = GOdata, classicFisher = resultFisher, 
                               topNodes = length(topGO::usedGO(GOdata)))
  goResults$classicFisher <- gsub("<", "", goResults$classicFisher)
  goResults$bonferroni <- stats::p.adjust(p = goResults$classicFisher, 
                                          method = "bonferroni")
  goResults$bh <- stats::p.adjust(goResults$classicFisher, 
                                  method = "BH")
  goResults$foldEnrichment <- round(goResults$Significant/goResults$Expected, 
                                    2)
  return(goResults)
}


runTopGO_comprehensive <- function(species, backgroundGenes, targetedGenes){
  bp_output <- runTopGO(ontology =  'BP', species = species, backgroundGenes =  backgroundGenes, targetedGenes = targetedGenes) %>% 
    mutate(ontology = "BP")
  mf_output <- runTopGO(ontology =  'MF', species = species, backgroundGenes =  backgroundGenes, targetedGenes = targetedGenes) %>% 
    mutate(ontology = "MF")
  cc_output <- runTopGO(ontology =  'CC', species = species, backgroundGenes =  backgroundGenes, targetedGenes = targetedGenes) %>% 
    mutate(ontology = "CC")
  combined_output <- bind_rows(bp_output, mf_output, cc_output)
  return(combined_output)
}

