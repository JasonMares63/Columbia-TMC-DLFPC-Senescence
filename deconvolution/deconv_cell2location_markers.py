import os
import glob
import torch
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from argparse import ArgumentParser

# https://github.com/BayraktarLab/cell2location
import cell2location
from cell2location.models import RegressionModel
from cell2location.utils.filtering import filter_genes

'''
For a given donor, run Cell2location in two steps:
- Infer reference gene expr. signatures for each cell type using only snRNA-seq data from donor.
- Perform cell type inference using the Cell2location model across all Visium arrays from said
  donor at once, using Visium array as batch index.
'''

parser = ArgumentParser()
parser.add_argument('-d', '--donor-id', type=str, required=True, 
	help='Donor ID (CUID) for snRNA-seq/Visium data')
parser.add_argument('-l', '--load-reference', action="store_true", default=False,
	help='Load snRNA-seq reference from pre-trained model in outputs directory -- otherwise, train model before deconvolution')
parser.add_argument('-f', '--full-data', action="store_true", default=False,
	help='Use reference cells from all donors in deconvolution (default to donor-matched cells).')
parser.add_argument('-n', '--num-cells', type=float, default=5,
	help='Mean number of cells per spot.')
args = parser.parse_args()

TRAIN_REF = not args.load_reference
USE_GPU = torch.cuda.is_available()
DONOR_MATCH = not args.full_data
N_CELLS = args.num_cells

data_dir = '/mnt/home/adaly/ceph/datasets/U54_BA9'
sc_dir = os.path.join(data_dir, 'snrna_anndata')
st_dir = os.path.join(data_dir, 'splotch_anndata')
if DONOR_MATCH:
	out_dir = os.path.join(data_dir, 'deconvolution', 'output_cell2location_markers')
else:
	out_dir = os.path.join(data_dir, 'deconvolution', 'output_cell2location_markers_fullref')
if not os.path.exists(out_dir):
	os.mkdir(out_dir)


# Map CU Donor ID to CGND Donor ID
subj_meta = pd.read_csv(os.path.join(sc_dir, 'snRNA-multiome-subjectlist-01092023.tsv'),
	header=0, index_col=0, sep='\t')
cu_id = args.donor_id
cgnd_id = subj_meta.loc[cu_id, 'Subject ID']
print('Donor: %s (%s)' % (cu_id, cgnd_id))

# Create a sub-directory for all deconvolution results from this donor
if not os.path.exists(os.path.join(out_dir, cu_id)):
	os.mkdir(os.path.join(out_dir, cu_id))


# Read snRNA-seq data
adata_sc = sc.read_h5ad(os.path.join(sc_dir, 'adata_U54_BA9_snrna_sennet_mapped_clusters.h5ad'))
adata_sc = adata_sc.raw.to_adata()
adata_sc.var.index = adata_sc.var.gene_name
celltype_label = 'mapped_cluster'


# Read ST data
adata_st = sc.read_h5ad(os.path.join(st_dir, 'adata_U54_BA9_splotch_lambdas.h5ad'))
donor_col = 'Level 3'

# Subset snRNA-seq/Visium data to donor being considered
if DONOR_MATCH:
	adata_sc = adata_sc[adata_sc.obs['cu_core_id']==cu_id, :].copy()
print(adata_sc.shape[0], 'nuclei')
adata_st = adata_st[adata_st.obs[donor_col]==cgnd_id].copy()
adata_st.obs['array'] = [x.split('/')[1] for x in adata_st.obs.index]
print(adata_st.shape[0], 'spots', len(adata_st.obs.array.unique()), 'arrays')

# Switch Visium indexing to common gene names to match snRNA-seq data
df_map = pd.read_csv(os.path.join(data_dir, 'gene_symbols.tsv'), header=0, index_col=0, sep='\t')
adata_st.var = adata_st.var.join(df_map, how='left')
adata_st = adata_st[:, adata_st.var.dropna().index].copy()
adata_st.var['ENSEMBL'] = adata_st.var.index
adata_st.var.index = adata_st.var.gene_name

# Load marker genes & subset both modalities
df_markers = pd.read_csv('U54_BA9_snrna_mapped_cluster_markers_balanced.csv', header=0, index_col=0)
genes_shared = df_markers['names'].unique()
adata_sc = adata_sc[:, genes_shared]
adata_st = adata_st[:, genes_shared]

# Load raw count data (instead of lambdas) in X
adata_st.X = adata_st.layers['X_counts'].todense()


# Perform reference signature inference
if DONOR_MATCH:
	ref_file = os.path.join(out_dir, cu_id, '%s_cell2loc_scref.h5ad' % cu_id)
else:
	ref_file = os.path.join(out_dir, 'cell2loc_scref.h5ad')

# Train a NB model on raw snRNA-seq data to estimate reference signatures (and save results)
if TRAIN_REF:
	cell2location.models.RegressionModel.setup_anndata(adata=adata_sc, labels_key=celltype_label)

	mod = RegressionModel(adata_sc)
	batch_size = 2500
	mod.train(max_epochs=250, batch_size=batch_size, use_gpu=USE_GPU)

	adata_sc = mod.export_posterior(
	    adata_sc,
	    sample_kwargs={'num_samples': 1000, 'batch_size': batch_size, 'use_gpu': USE_GPU}
	)
	adata_sc.write(ref_file)
	mod.plot_history(20)
	mod.plot_QC()
	plt.close()
# Load previously saved model
else:
	adata_sc = sc.read_h5ad(ref_file)

# Export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_sc.varm.keys():
	print('yay')
	inf_aver = adata_sc.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' \
	for i in adata_sc.uns['mod']['factor_names']]].copy()
else:
	inf_aver = adata_sc.var[[f'means_per_cluster_mu_fg_{i}' \
	for i in adata_sc.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_sc.uns['mod']['factor_names']


# Perform deconvolution of Visium data
cell2location.models.Cell2location.setup_anndata(adata=adata_st, batch_key="array")

mod = cell2location.models.Cell2location(
	adata_st, cell_state_df=inf_aver,
	N_cells_per_location=N_CELLS,
	# hyperparameter controlling normalisation of
	# within-experiment variation in RNA detection:
	detection_alpha=20
)
mod.train(max_epochs=30000,
	# train using full data (batch_size=None)
	batch_size=None,
	# use all data points in training because
	# we need to estimate cell abundance at all locations
	train_size=1,
	use_gpu=USE_GPU)
mod.plot_history(1000)
plt.close()

# Export deconvolution results
adata_st = mod.export_posterior(
    adata_st, 
    sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': USE_GPU}
)
deconv_file = os.path.join(out_dir, cu_id, '%s_cell2loc_deconv.h5ad' % cu_id)
adata_st.write(deconv_file)

fig = mod.plot_QC()
plt.savefig(os.path.join(out_dir, cu_id, '%s_cell2loc_deconv_reconst.png' % cu_id))
plt.close()

fig = mod.plot_spatial_QC_across_batches()
plt.savefig(os.path.join(out_dir, cu_id, '%s_cell2loc_deconv_batch.png' % cu_id))
plt.close()

