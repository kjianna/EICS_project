###############################################################
## EICS Definition
###############################################################

library(Seurat)
library(here)
library(dplyr)
library(ggplot2)

###############################################################
## Load data
###############################################################

fetal_EC <- readRDS(
  here(
    "results",
    "reference",
    "fetal_EC_module_score.rds"
  )
)

###############################################################
## Build EICS matrix
###############################################################

EICS_matrix <-
  fetal_EC@meta.data %>%
  dplyr::select(
    Paracrine_score,
    Mechanical_score,
    Metabolic_score
  )

head(EICS_matrix)

###############################################################
## PCA on Instruction Space
###############################################################

EICS_pca <-
  prcomp(
    EICS_matrix,
    center = TRUE,
    scale. = TRUE
  )

EICS_embedding <-
  as.data.frame(EICS_pca$x)

head(EICS_embedding)

EICS_embedding$cluster <-
  fetal_EC$seurat_clusters

EICS_embedding$week <-
  fetal_EC$week

EICS_embedding$region <-
  fetal_EC$region

ggplot(
  EICS_embedding,
  aes(
    PC1,
    PC2,
    color = cluster
  )
) +
  geom_point(size = 2) +
  theme_classic() +
  labs(
    title = "Instruction Space",
    x = "PC1",
    y = "PC2"
  )

fetal_EC$week <- factor(
  fetal_EC$week,
  levels = c(
    "5W","6W","7W","9W",
    "10W","13W","15W",
    "17W","20W","22W",
    "23W","24W","25W"
  )
)

ggplot(
  EICS_embedding,
  aes(
    PC1,
    PC2,
    color = week
  )
) +
  geom_point(size = 2) +
  theme_classic()

ggplot(
  EICS_embedding,
  aes(
    PC1,
    PC2,
    color = region
  )
) +
  geom_point(size = 2) +
  theme_classic()

###############################################################
## PCA loading
###############################################################

EICS_pca$rotation
summary(EICS_pca)

library(plotly)

plot_ly(
  EICS_matrix,
  x = ~Paracrine_score,
  y = ~Mechanical_score,
  z = ~Metabolic_score,
  color = fetal_EC$seurat_clusters,
  colors = "Set1",
  type = "scatter3d",
  mode = "markers"
)