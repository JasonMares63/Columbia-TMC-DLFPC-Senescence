library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(patchwork)

# See "save_anndata_rds.R" for how to convert AnnData to Seurat RDS
ba9 = LoadSeuratRds("adata_ba9_visium_rawcounts.rds")

# 1. Split dataset by array; log-transform & select HVGs separately

split.by = "visarray"
# use all arrays from one individual (SD019/15, Male (50yo)) as a reference
reference.names = c("CU001-U54-HRA-010-C", "CU001-U54-HRA-010-D",
                    "CU002-U54-HRA-010-C", "CU002-U54-HRA-010-D")

# 2. Select integration features across all arrays

ba9.list = SplitObject(ba9, split.by=split.by)
ba9.list = lapply(X=ba9.list, FUN=function(x){
  x = NormalizeData(x)
  x = FindVariableFeatures(x, selection.method="vst", nfeatures=2000)
})

# 3. Scale integration features and run PCA separately on each array

features = SelectIntegrationFeatures(object.list=ba9.list)
ba9.list = lapply(X=ba9.list, FUN = function(x){
  x = ScaleData(x, features=features, verbose=FALSE)
  x = RunPCA(x, features=features, verbose=FALSE)
})

# 4. Use RPCA to find integration anchors with all arrays from reference individual

reference = match(reference.names, names(ba9.list))
anchors = FindIntegrationAnchors(ba9.list, reference=reference,
                                 reduction="rpca", dims=1:50)
ba9.integrated = IntegrateData(anchorset=anchors, dims=1:50)

## 5. Integrate data, scale, then run PCA and UMAP.

#ba9.integrated = ScaleData(ba9.integrated, verbose=FALSE)
#ba9.integrated = RunPCA(ba9.integrated, verbose=FALSE)
#ba9.integrated = RunUMAP(ba9.integrated, dims=1:50)

saveRDS(ba9.integrated, file="ba9_array_integrated.rds")