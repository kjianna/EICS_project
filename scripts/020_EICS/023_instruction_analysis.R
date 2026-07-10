library(Seurat)
library(here)
library(ggplot2)

fetal_EC <- readRDS(
  here(
    "results",
    "reference",
    "fetal_EC_module_score.rds"
  )
)

VlnPlot(
  fetal_EC,
  features = c(
    "Paracrine_score",
    "Mechanical_score",
    "Metabolic_score"
  ),
  pt.size = 0
)

grep(
  "Para|Mech|Meta",
  colnames(fetal_EC@meta.data),
  value = TRUE
)

library(dplyr)

cluster_score <-
  fetal_EC@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    Paracrine = mean(Paracrine_score),
    Mechanical = mean(Mechanical_score),
    Metabolic = mean(Metabolic_score)
  )

kruskal.test(
  Paracrine_score ~ seurat_clusters,
  data = fetal_EC@meta.data
)

kruskal.test(
  Mechanical_score ~ seurat_clusters,
  data = fetal_EC@meta.data
)

kruskal.test(
  Metabolic_score ~ seurat_clusters,
  data = fetal_EC@meta.data
)

pairwise.wilcox.test(
  fetal_EC$Paracrine_score,
  fetal_EC$seurat_clusters,
  p.adjust.method = "BH"
)

pairwise.wilcox.test(
  fetal_EC$Mechanical_score,
  fetal_EC$seurat_clusters,
  p.adjust.method = "BH"
)

pairwise.wilcox.test(
  fetal_EC$Metabolic_score,
  fetal_EC$seurat_clusters,
  p.adjust.method = "BH"
)

cluster_score

cluster_summary <-
  fetal_EC@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    Paracrine_mean = mean(Paracrine_score),
    Paracrine_sd = sd(Paracrine_score),
    
    Mechanical_mean = mean(Mechanical_score),
    Mechanical_sd = sd(Mechanical_score),
    
    Metabolic_mean = mean(Metabolic_score),
    Metabolic_sd = sd(Metabolic_score)
  )

cluster_summary

fetal_EC@meta.data %>%
  count(seurat_clusters)

library(pheatmap)

score_matrix <-
  as.matrix(cluster_score[,2:4])

rownames(score_matrix) <-
  paste0("Cluster_", cluster_score$seurat_clusters)

pheatmap(score_matrix)

pheatmap(
  score_matrix,
  scale = "column",
  main = "Scaled Instruction Scores"
)

###############################################################
## Week comparison
###############################################################

fetal_EC$week <- factor(
  fetal_EC$week,
  levels = c(
    "5W","6W","7W","9W",
    "10W","13W","15W",
    "17W","20W","22W",
    "23W","24W","25W"
  )
)

VlnPlot(
  fetal_EC,
  features = "Paracrine_score",
  group.by = "week",
  pt.size = 0
)

VlnPlot(
  fetal_EC,
  features = "Mechanical_score",
  group.by = "week",
  pt.size = 0
)

VlnPlot(
  fetal_EC,
  features = "Metabolic_score",
  group.by = "week",
  pt.size = 0
)

library(dplyr)

week_summary <-
  fetal_EC@meta.data %>%
  group_by(week) %>%
  summarise(
    Paracrine_mean = mean(Paracrine_score),
    Mechanical_mean = mean(Mechanical_score),
    Metabolic_mean = mean(Metabolic_score),
    
    Paracrine_sd = sd(Paracrine_score),
    Mechanical_sd = sd(Mechanical_score),
    Metabolic_sd = sd(Metabolic_score),
    
    n = n()
  )

week_summary

library(tidyr)
library(ggplot2)

week_plot <-
  week_summary %>%
  pivot_longer(
    cols = c(
      Paracrine_mean,
      Mechanical_mean,
      Metabolic_mean
    ),
    names_to = "Instruction",
    values_to = "Score"
  )

week_plot$Instruction <-
  factor(
    week_plot$Instruction,
    levels = c(
      "Paracrine_mean",
      "Mechanical_mean",
      "Metabolic_mean"
    ),
    labels = c(
      "Paracrine",
      "Mechanical",
      "Metabolic"
    )
  )

ggplot(
  week_plot,
  aes(
    x = week,
    y = Score,
    color = Instruction,
    group = Instruction
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_classic()

kruskal.test(
  Paracrine_score ~ week,
  data = fetal_EC@meta.data
)

kruskal.test(
  Mechanical_score ~ week,
  data = fetal_EC@meta.data
)

kruskal.test(
  Metabolic_score ~ week,
  data = fetal_EC@meta.data
)

colnames(fetal_EC@meta.data)

###############################################################
## Sample composition
###############################################################

fetal_EC@meta.data %>%
  count(sample)

table(fetal_EC$sample, fetal_EC$week)

table(fetal_EC$region)
levels(as.factor(fetal_EC$region))

###############################################################
## Major cardiac regions only
###############################################################

major_region <- c(
  "RA",
  "LA",
  "RV",
  "LV"
)

fetal_EC_region <-
  subset(
    fetal_EC,
    subset = region %in% major_region
  )

VlnPlot(
  fetal_EC_region,
  features = "Paracrine_score",
  group.by = "region",
  pt.size = 0
)

VlnPlot(
  fetal_EC_region,
  features = "Mechanical_score",
  group.by = "region",
  pt.size = 0
)

VlnPlot(
  fetal_EC_region,
  features = "Metabolic_score",
  group.by = "region",
  pt.size = 0
)

region_summary <-
  fetal_EC_region@meta.data %>%
  group_by(region) %>%
  summarise(
    Paracrine_mean = mean(Paracrine_score),
    Mechanical_mean = mean(Mechanical_score),
    Metabolic_mean = mean(Metabolic_score),
    
    Paracrine_sd = sd(Paracrine_score),
    Mechanical_sd = sd(Mechanical_score),
    Metabolic_sd = sd(Metabolic_score),
    
    n = n()
  )

region_summary

kruskal.test(
  Paracrine_score ~ region,
  data = fetal_EC_region@meta.data
)

kruskal.test(
  Mechanical_score ~ region,
  data = fetal_EC_region@meta.data
)

kruskal.test(
  Metabolic_score ~ region,
  data = fetal_EC_region@meta.data
)

###############################################################
## Correlation analysis
###############################################################

cor(
  fetal_EC$Paracrine_score,
  fetal_EC$Mechanical_score
)

cor(
  fetal_EC$Paracrine_score,
  fetal_EC$Metabolic_score
)

cor(
  fetal_EC$Mechanical_score,
  fetal_EC$Metabolic_score
)

FeatureScatter(
  fetal_EC,
  feature1 = "Paracrine_score",
  feature2 = "Mechanical_score"
)

FeatureScatter(
  fetal_EC,
  feature1="Paracrine_score",
  feature2="Metabolic_score"
)

FeatureScatter(
  fetal_EC,
  feature1="Mechanical_score",
  feature2="Metabolic_score"
)