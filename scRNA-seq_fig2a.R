#identify 8CLCs in naive hESC cultures

library(Seurat)
library(dplyr)
library(ggplot2)
library(GEOquery)

#load data
geo <- getGEO("GSE232939", GSEMatrix =TRUE)
exprSet <- exprs(geo$GSE232939_series_matrix.txt.gz)
meta <- pData(geo$GSE232939_series_matrix.txt.gz)
meta <- meta %>% select(title, `cell type:ch1`, `culture condition:ch1`)
colnames(meta) <- c("cell", "cell_type", "culture_condition")
#prepare seurat object
seurat_obj <- CreateSeuratObject(counts = exprSet, meta.data = meta)