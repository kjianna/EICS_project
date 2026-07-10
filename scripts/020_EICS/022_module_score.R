############################################################
## Project : EICS
## Script  : 02_module_score.R
## Purpose :
##   Calculate instruction module scores
##   using literature-derived gene sets.
############################################################

library(Seurat)
library(tidyverse)
library(readr)
library(here)

############################################################
## Load Seurat object
############################################################

'''
fetal_EC <- readRDS(
  here(
    "results",
    "reference",
    "fetal_EC_reclustered.rds"
  )
)
'''

dir.create(
  here("results", "reference"),
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  fetal_EC,
  here(
    "results",
    "reference",
    "fetal_EC_module_score.rds"
  )
)

############################################################
## Load gene sets
############################################################

paracrine <- read_csv(
  here(
    "gene_set",
    "EICS_v1",
    "paracrine_instruction.csv"
  ),
  show_col_types = FALSE
)

mechanical <- read_csv(
  here(
    "gene_set",
    "EICS_v1",
    "mechanical_instruction.csv"
  ),
  show_col_types = FALSE
)

metabolic <- read_csv(
  here(
    "gene_set",
    "EICS_v1",
    "metabolic_instruction.csv"
  ),
  show_col_types = FALSE
)

############################################################
## Filter genes
############################################################

paracrine_genes <- intersect(
  paracrine$Gene,
  rownames(fetal_EC)
)

mechanical_genes <- intersect(
  mechanical$Gene,
  rownames(fetal_EC)
)

metabolic_genes <- intersect(
  metabolic$Gene,
  rownames(fetal_EC)
)

length(paracrine_genes)
length(mechanical_genes)
length(metabolic_genes)

############################################################
## Calculate module scores
############################################################

fetal_EC <- AddModuleScore(
  fetal_EC,
  features = list(paracrine_genes),
  name = "Paracrine"
)

fetal_EC <- AddModuleScore(
  fetal_EC,
  features = list(mechanical_genes),
  name = "Mechanical"
)

fetal_EC <- AddModuleScore(
  fetal_EC,
  features = list(metabolic_genes),
  name = "Metabolic"
)

colnames(paracrine)
head(paracrine)

############################################################
## Check scores
############################################################

head(
  fetal_EC@meta.data %>%
    select(
      starts_with("Paracrine"),
      starts_with("Mechanical"),
      starts_with("Metabolic")
    )
)

############################################################
## Save object
############################################################

saveRDS(
  fetal_EC,
  here(
    "results",
    "reference",
    "fetal_EC_module_score.rds"
  )
)

############################################################
## Rename score columns
############################################################

fetal_EC$Paracrine_score <- fetal_EC$Paracrine1
fetal_EC$Mechanical_score <- fetal_EC$Mechanical1
fetal_EC$Metabolic_score <- fetal_EC$Metabolic1