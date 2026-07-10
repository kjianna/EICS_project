############################################################
## Project : EICS
## Script  : 01_instruction_panel.R
## Purpose :
##   Build literature-derived instruction gene sets
############################################################

library(tidyverse)
library(here)

############################################################
## Paracrine instruction genes
############################################################

paracrine <- tibble(
  Gene = c(
    "NRG1","IGF2",
    "BMP2","BMP4","BMP10",
    "WNT2","WNT5A",
    "FGF9","FGF16",
    "INHBA","INHBB",
    "DLL4","JAG1","JAG2",
    "EFNB2","EFNA1",
    "CXCL12",
    "PDGFB",
    "EDN1",
    "HBEGF",
    "VEGFA"
  ),
  Category = "Paracrine",
  Priority = c(
    "Core","Core",
    "Core","Core","Core",
    "Core","Support",
    "Core","Support",
    "Candidate","Support",
    "Core","Core","Support",
    "Core","Support",
    "Core",
    "Core",
    "Core",
    "Support",
    "Core"
  )
)

############################################################
## Mechanical instruction genes
############################################################

mechanical <- tibble(
  Gene = c(
    "KLF2","KLF4",
    "PIEZO1",
    "NOS3",
    "YAP1","WWTR1","TEAD1",
    "GJA1",
    "ENG"
  ),
  Category = "Mechanical",
  Priority = c(
    "Core","Support",
    "Core",
    "Core",
    "Core","Core","Core",
    "Support",
    "Support"
  )
)

############################################################
## Metabolic instruction genes
############################################################

metabolic <- tibble(
  Gene = c(
    "FABP4",
    "CD36",
    "CA4",
    "GPIHBP1",
    "PLVAP",
    "ACE",
    "APOD",
    "IGFBP7",
    "RAMP2",
    "ADM"
  ),
  Category = "Metabolic",
  Priority = c(
    "Core","Core","Core","Core","Core",
    "Core","Support","Support","Support","Support"
  )
)

############################################################
## Save gene sets
############################################################

write_csv(
  paracrine,
  here("gene_sets","EICS_v1","paracrine_instruction.csv")
)

write_csv(
  mechanical,
  here("gene_sets","EICS_v1","mechanical_instruction.csv")
)

write_csv(
  metabolic,
  here("gene_sets","EICS_v1","metabolic_instruction.csv")
)

message("EICS v1 gene sets saved.")############################################################
## Project : EICS
## Script  : 01_instruction_panel.R
## Purpose :
##   Build literature-derived instruction gene sets
############################################################

library(tidyverse)
library(here)

############################################################
## Paracrine instruction genes
############################################################

paracrine <- tibble(
  Gene = c(
    "NRG1","IGF2",
    "BMP2","BMP4","BMP10",
    "WNT2","WNT5A",
    "FGF9","FGF16",
    "INHBA","INHBB",
    "DLL4","JAG1","JAG2",
    "EFNB2","EFNA1",
    "CXCL12",
    "PDGFB",
    "EDN1",
    "HBEGF",
    "VEGFA"
  ),
  Category = "Paracrine",
  Priority = c(
    "Core","Core",
    "Core","Core","Core",
    "Core","Support",
    "Core","Support",
    "Candidate","Support",
    "Core","Core","Support",
    "Core","Support",
    "Core",
    "Core",
    "Core",
    "Support",
    "Core"
  )
)

############################################################
## Mechanical instruction genes
############################################################

mechanical <- tibble(
  Gene = c(
    "KLF2","KLF4",
    "PIEZO1",
    "NOS3",
    "YAP1","WWTR1","TEAD1",
    "GJA1",
    "ENG"
  ),
  Category = "Mechanical",
  Priority = c(
    "Core","Support",
    "Core",
    "Core",
    "Core","Core","Core",
    "Support",
    "Support"
  )
)

############################################################
## Metabolic instruction genes
############################################################

metabolic <- tibble(
  Gene = c(
    "FABP4",
    "CD36",
    "CA4",
    "GPIHBP1",
    "PLVAP",
    "ACE",
    "APOD",
    "IGFBP7",
    "RAMP2",
    "ADM"
  ),
  Category = "Metabolic",
  Priority = c(
    "Core","Core","Core","Core","Core",
    "Core","Support","Support","Support","Support"
  )
)

############################################################
## Save gene sets
############################################################

write_csv(
  paracrine,
  here("gene_sets","EICS_v1","paracrine_instruction.csv")
)

write_csv(
  mechanical,
  here("gene_sets","EICS_v1","mechanical_instruction.csv")
)

write_csv(
  metabolic,
  here("gene_sets","EICS_v1","metabolic_instruction.csv")
)

message("EICS v1 gene sets saved.")