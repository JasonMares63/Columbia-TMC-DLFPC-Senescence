library(Seurat)
library(ggplot2)
library(patchwork)

ba9.integrated = readRDS("ba9_array_integrated.rds")

# Read in tissue annotation info (BA46 vs BA9)
df.meta = read.csv("../U54_Metadata_Full.tsv", header=1, sep="\t")
res = sapply(df.meta$tissue_id, FUN=function(x){
  parts = strsplit(x, '_')
  first_elements = sapply(parts, FUN=function(x) x[1])
  return(first_elements)
})
df.meta$tissue_type = res

# Join sample metadata to spot metadata from ba9.integrated
spot.meta = ba9.integrated@meta.data
res = merge(x=spot.meta, y=df.meta, by.x="visarray", by.y="library_sample_id", 
            all.x=TRUE)
ba9.integrated = AddMetaData(object=ba9.integrated, metadata=res$tissue_type, 
            col.name="tissue.type")

df.umap = Embeddings(ba9.integrated, reduction="umap")
df.pca = Embeddings(ba9.integrated, reduction="pca")[,1:3]

write.csv(df.umap, file="ba9_array_integrated_umap.csv")

aar_colors = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2")
DimPlot(ba9.integrated, group.by="region", raster=FALSE, pt.size=0.001, alpha=0.7, cols=aar_colors)
ggsave(path="plots", filename="ba9_umap_region.png", device="png",
       units="in", width=12, height=10)

DimPlot(ba9.integrated, group.by="Level.1", raster=FALSE)
ggsave(path="plots", filename="ba9_umap_l1.png", device="png",
       units="in", width=12, height=10)

DimPlot(ba9.integrated, group.by="Level.2", raster=FALSE)
ggsave(path="plots", filename="ba9_umap_l2.png", device="png",
       units="in", width=12, height=10)

DimPlot(ba9.integrated, group.by="Level.3", raster=FALSE)
ggsave(path="plots", filename="ba9_umap_l3.png", device="png",
       units="in", width=13, height=10)

DimPlot(ba9.integrated, group.by="tissue.type", raster=FALSE)
ggsave(path="plots", filename="ba9_umap_tissue.png", device="png",
       units="in", width=12, height=10)