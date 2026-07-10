####################################################
## Milestone 1. Create Seurat object
####################################################

library(Seurat)
library(data.table)
library(readxl)
library(dplyr)
library(here)

sessionInfo()

# expression matrix

counts <- read.table(
  here(
    "data",
    "raw",
    "fetal",
    "Cui_2019_GSE106118",
    "expression",
    "GSE106118_UMI_count_merge.txt.gz"
  ),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

dim(counts)

head(counts[,1:5])

gene_names <- counts$Gene

counts <- counts[, -1]

rownames(counts) <- gene_names

counts <- as.matrix(counts)

dim(counts)

counts[1:5,1:5]

rownames(counts)[1:10]

colnames(counts)[1:10]

# Seurat object

fetal <- CreateSeuratObject(
  counts = counts,
  project = "Cui2019",
  min.cells = 3,
  min.features = 200
)

fetal
dim(fetal)
head(fetal@meta.data)

####################################################
## Milestone 2. Add cell annotation
####################################################

#annotation

annotation <- readxl::read_excel(
  here(
    "data",
    "raw",
    "fetal",
    "Cui_2019_GSE106118",
    "metadata",
    "Table_S2_cell_annotation.xlsx"
  )
)

dim(annotation)
head(annotation)
colnames(annotation)

colnames(annotation) <- c("cell", "cluster")
head(annotation)
head(colnames(fetal))

library(tidyr)

annotation <- annotation %>%
  mutate(
    cluster_id = sub(" .*", "", cluster),
    cell_type = gsub(".*\\((.*)\\)", "\\1", cluster)
  )

cluster_map <- c(
  "C1" = "Early_Heart",
  "C2" = "Cardiomyocyte",
  "C3" = "Fibroblast",
  "C4" = "Endothelial",
  "C5" = "Valvular",
  "C6" = "Epicardium",
  "C7" = "Mast",
  "C8" = "Macrophage",
  "C9" = "T_B_cell"
)

annotation <- annotation %>%
  mutate(
    cell_type = cluster_map[cluster_id]
  )

head(annotation)

# Seurat metadata

annotation_df <- as.data.frame(annotation)

rownames(annotation_df) <- annotation_df$cell
annotation_df$cell <- NULL

head(annotation_df)

# Seurat object + metadata

fetal <- AddMetaData(
  object = fetal,
  metadata = annotation_df
)

head(fetal@meta.data)

table(fetal$cluster_id)

table(fetal$cell_type)

sum(!is.na(fetal$cell_type))

fetal_ref <- subset(
  fetal,
  cells = rownames(annotation_df)
)

fetal_ref

table(fetal_ref$cell_type)

####################################################
## Milestone 3. Build endothelial reference
####################################################

# Extract EC

fetal_EC <- subset(
  fetal_ref,
  subset = cell_type == "Endothelial"
)

fetal_EC
table(fetal_EC$cell_type)

ec_subtype <- readxl::read_excel(
  here(
    "data",
    "raw",
    "fetal",
    "Cui_2019_GSE106118",
    "metadata",
    "Table_S5_EC_subtype.xlsx"
  ),
  sheet = 2
)

dim(ec_subtype)
head(ec_subtype)
colnames(ec_subtype)

fetal_EC <- subset(
  fetal_ref,
  subset = cell_type == "Endothelial"
)

fetal_EC

sum(colnames(fetal_EC) %in% ec_subtype$Cells)

length(intersect(
  colnames(fetal_EC),
  ec_subtype$Cells
))

setdiff(
  ec_subtype$Cells,
  colnames(fetal_EC)
)[1:20]

############################################################
## STEP 4. Reprocess endothelial cells
############################################################

fetal_EC <- NormalizeData(fetal_EC)

fetal_EC <- FindVariableFeatures(
  fetal_EC,
  selection.method = "vst",
  nfeatures = 2000
)

fetal_EC <- ScaleData(fetal_EC)

fetal_EC <- RunPCA(
  fetal_EC,
  npcs = 30
)

ElbowPlot(fetal_EC)
print(fetal_EC)
DimPlot(fetal_EC, reduction = "pca")

############################################################
## STEP 5. Clustering endothelial cells
############################################################

fetal_EC <- FindNeighbors(
  fetal_EC,
  dims = 1:10
)

fetal_EC <- FindClusters(
  fetal_EC,
  resolution = 0.3
)

fetal_EC <- RunUMAP(
  fetal_EC,
  dims = 1:10
)

DimPlot(
  fetal_EC,
  reduction = "umap",
  label = TRUE
)

table(Idents(fetal_EC))

cell_info <- data.frame(
  cell = colnames(fetal_ref)
)

parts <- strsplit(cell_info$cell, "_")

cell_info$week <- sapply(parts, function(x) gsub("HE", "", x[1]))
cell_info$sample <- sapply(parts, function(x) x[2])
cell_info$region <- sapply(parts, function(x) sub("\\..*", "", x[3]))

head(cell_info)

rownames(cell_info) <- cell_info$cell
cell_info$cell <- NULL

fetal_ref <- AddMetaData(
  fetal_ref,
  metadata = cell_info
)

fetal_EC <- subset(
  fetal_ref,
  subset = cell_type == "Endothelial"
)

head(fetal_EC@meta.data)

fetal_EC <- NormalizeData(fetal_EC)

fetal_EC <- FindVariableFeatures(
  fetal_EC,
  selection.method = "vst",
  nfeatures = 2000
)

fetal_EC <- ScaleData(fetal_EC)

fetal_EC <- RunPCA(fetal_EC)

fetal_EC <- FindNeighbors(
  fetal_EC,
  dims = 1:10
)

fetal_EC <- FindClusters(
  fetal_EC,
  resolution = 0.3
)

fetal_EC <- RunUMAP(
  fetal_EC,
  dims = 1:10
)

DimPlot(
  fetal_EC,
  reduction = "umap",
  group.by = "week"
)

DimPlot(
  fetal_EC,
  group.by = "region"
)

DimPlot(
  fetal_EC,
  group.by = "sample"
)

markers <- FindAllMarkers(
  fetal_EC,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

library(dplyr)

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

top10

print(top10, n = 40)

write.csv(
  markers,
  "results/fetal_EC_markers.csv",
  row.names = FALSE
)

View(top10)

FeaturePlot(
  fetal_EC,
  features = c(
    "PECAM1",
    "VWF",
    "KDR",
    "FLT1",
    "EMCN",
    "CDH5"
  ),
  reduction = "umap",
  ncol = 3
)

FeaturePlot(
  fetal_EC,
  features = c(
    "SOX17",
    "EFNB2",
    "GJA5",
    "NR2F2",
    "NPR3",
    "PLVAP"
  ),
  reduction = "umap",
  ncol = 3
)

markers <- FindAllMarkers(
  fetal_EC,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  markers,
  "results/fetal_EC_all_markers.csv",
  row.names = FALSE
)

cluster0 <- markers %>%
  filter(cluster == 0) %>%
  arrange(desc(avg_log2FC))

cluster1 <- markers %>%
  filter(cluster == 1) %>%
  arrange(desc(avg_log2FC))

cluster2 <- markers %>%
  filter(cluster == 2) %>%
  arrange(desc(avg_log2FC))

cluster3 <- markers %>%
  filter(cluster == 3) %>%
  arrange(desc(avg_log2FC))