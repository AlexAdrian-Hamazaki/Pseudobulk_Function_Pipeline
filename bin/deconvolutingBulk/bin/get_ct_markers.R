library(Seurat)
#devtools::install_github('satijalab/seurat-data', force = TRUE)
library(SeuratData)
library(ggplot2)
#InstallData('ifnb')
LoadData('ifnb')

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

# Normalize the data
ifnb <- NormalizeData(ifnb)

# Find DE features between CD16 Mono and CD1 Mono
Idents(ifnb) <- "seurat_annotations"

head(Idents(ifnb))

monocyte.de.markers <- FindMarkers(ifnb, ident.1 = "CD16 Mono", ident.2 = "CD14 Mono")

head(monocyte.de.markers)

