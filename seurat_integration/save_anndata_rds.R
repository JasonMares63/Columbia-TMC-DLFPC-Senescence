library(Matrix)
library(Seurat)
library(SeuratObject)
library(reticulate)  # allows use of Python functions in R; for some reason having trouble with Rstudio using it

sc = import("scanpy")

# For some reason, when reticulate reads the CSR-formatted integer matrix in adata.X,
# it loads it into a dgRMatrix object, which should only hold "double" and throws errors on any subsequent operation.
#
# To work around, I saved the CSR matrix into .mtx format using scipy and read it in separately:
"
import scanpy as sc 
from scipy.io import mmwrite

adata = sc.read_h5ad('../adata_U54_BA9_visium_rawcounts.h5ad')
mmwrite('../adata_U54_BA9_visium_rawcounts.mtx', adata.X.tocsr())
"
data <- sc$read_h5ad("../adata_U54_BA9_visium_rawcounts.h5ad")
counts_matrix = readMM('../adata_U54_BA9_visium_rawcounts.mtx')

counts <- t(counts_matrix)
colnames(counts) <-  data$obs_names$to_list()
rownames(counts) <-  data$var_names$to_list()
counts <- Matrix::Matrix(as.matrix(counts), sparse = T)

seurat <- CreateSeuratObject(counts)
seurat <- AddMetaData(seurat,  data$obs)
SaveSeuratRds(seurat, 'adata_ba9_visium_rawcounts.rds')