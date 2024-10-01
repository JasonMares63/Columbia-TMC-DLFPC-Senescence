import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

adata = sc.read_h5ad('adata_U54_BA9_splotch_lambdas.h5ad')
print(adata, flush=True)

# Cap lambdas for each gene to their 99th percentile to remove effect of outliers                                                                                                 
cap_val = np.percentile(adata.X, 99, axis=0)
print('99th percentile:', cap_val, flush=True)

lambdas_clipped = np.clip(adata.X, np.zeros_like(cap_val), cap_val)

# Scale lambdas for each gene across all spots w/in each individual (account for biological variation)                                                                           
adata_scaled = ad.AnnData(lambdas_clipped, obs=adata.obs, var=adata.var, dtype=np.float32)

adatas_list = []

# Create new AnnData by separating spots by individual, scaling, then concatenating (to remove biological variation).                                                             
# NOTE: this will change the order of observations!                                                                                                                               
for indiv in adata.obs['Level 3'].unique():
    adata_indiv = adata_scaled[adata_scaled.obs['Level 3'] == indiv].copy()
    sc.pp.scale(adata_indiv, max_value=10.0)

    adatas_list.append(adata_indiv)

adata_scaled = ad.concat(adatas_list)
print(adata_scaled, flush=True)

# Ensure ordering of observations is consistent between original (raw) and scaled AnnDatas.                                                                                       
adata.layers['X_scaled'] = adata_scaled[adata.obs.index,:].X
print(adata, flush=True)

assert not np.isnan(adata.X).any(), 'Mean lambda values (adata.X) contain NaNs!'
assert not np.isnan(adata.layers['X_scaled']).any(), 'Scaled & clipped mean lambda values (X_scaled layer) contain NaNs!'

adata.write('adata_U54_BA9_splotch_lambdas_scaled.h5ad')
