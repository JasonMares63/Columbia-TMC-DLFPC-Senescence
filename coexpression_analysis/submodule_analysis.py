import os
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster


# Compute mean expression within a cell group
def grouped_obs_mean(adata, group_key, layer=None, gene_symbols=None):
	if layer is not None:
		getX = lambda x: x.layers[layer]
	else:
		getX = lambda x: x.X
	if gene_symbols is not None:
		new_idx = adata.var[gene_symbols].values
	else:
		new_idx = adata.var_names

	grouped = adata.obs.groupby(group_key, observed=False)
	
	out = {}
	for group, idx in grouped.indices.items():
		X = getX(adata[idx])
		out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))    
	out = pd.DataFrame(out, index=new_idx)
	return out

# Preprocess raw snRNA-seq count data for clustering
def preprocess_snrna(adata):
	adata_pp = sc.pp.normalize_total(adata, 1e6, copy=True)    # TPM normalization

	# discard lowly expressed genes (<10 TPM max across dataset)
	adata_pp = adata_pp[:, adata_pp.X.max(axis=0) > 10].copy()
	return adata_pp

# For a given geneset, find subclusters associated with given cell type annotations in snRNA-seq data
def plot_submodules(adata, geneset, obs_ctype, vmax=None, min_genes=10, display_genes=None):
	df_means = grouped_obs_mean(adata[:, geneset], obs_ctype)

	# Normalize gene expression between 0 and 1 across cell types
	df_means = df_means / df_means.max(axis=1).values[:, None] 

	# Determine cluster membership by distance cutoff
	Z = linkage(df_means.values, method='average', metric='cosine')
	max_d = 0.54 * max(Z[:,2])
	clusters = fcluster(Z, max_d, criterion='distance')

	# Map all clusters with >min_genes to a Tab20 color; all others to white.
	k = len(np.unique(clusters))
	cmap = plt.get_cmap('tab20')
	cNorm = colors.Normalize(vmin=0, vmax=20)
	scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmap)
	k2col = [scalarMap.to_rgba(np.remainder(i, 20)) if np.sum(clusters==i+1) > min_genes else (1,1,1,1) for i in range(k)]
    
	figsize=(10,10)
	bottom = 0.25

	# Render mean expression as a heatmap with rows grouped (and colored) by subcluster
	if vmax is None:
		vmax = df_means.values.max()
	cg = sns.clustermap(df_means, method='average', metric='cosine',
				   row_cluster=True, row_linkage=Z, row_colors=[k2col[i-1] for i in clusters],
				   col_cluster=False,
				   vmin=0, vmax=vmax,
				   figsize=figsize,
				   cmap='Greys',
				   cbar_kws={'orientation':'horizontal'})
	cg.ax_heatmap.set_xticks(np.arange(len(df_means.columns))+0.5)
	cg.ax_heatmap.set_xticklabels(df_means.columns, fontsize = 16)
	cg.figure.subplots_adjust(bottom=bottom)

    # Reduce number of genes shown so we can increase font size
	if display_genes is not None:
		row_inds_ordered = cg.dendrogram_row.reordered_ind
		row_genes_ordered = df_means.index[row_inds_ordered]
		yticks = [list(row_genes_ordered).index(g)+0.5 for g in display_genes]
		sort_inds = np.argsort(yticks)
		yticks = np.array(yticks)[sort_inds]
		yticklabels = np.array(display_genes)[sort_inds]
		ytl = []
		for i, (y,g) in enumerate(zip(yticks, yticklabels)):
			ytl.append(g)
			if i > 0 and np.abs(yticks[i-1]-y) < 30:
				ytl[i] += ', ' + ytl[i-1]
				ytl[i-1] = ''
		cg.ax_heatmap.set_yticks(yticks)
		cg.ax_heatmap.set_yticklabels(ytl, fontsize=14)
		print(yticks, yticklabels)
	else:
		cg.ax_heatmap.set_yticks([])
	cg.ax_heatmap.set_ylabel('')
    
	x0, _y0, _w, _h = cg.cbar_pos
	cg.ax_cbar.set_position([x0, 0.9, cg.ax_row_dendrogram.get_position().width+0.05, 0.02])
	cg.ax_cbar.set_title('Scaled mean expression', fontsize=14)
	cg.ax_cbar.tick_params(axis='x', length=10)
	for t in cg.ax_cbar.get_xticklabels():
		t.set_fontsize(14)

	df_label = pd.DataFrame({'submodule': clusters}, index=df_means.index)

	# Re-label all submodules s.t. only those with >min_genes receive an index; all others are NaN
	smod_cnt = 0
	for lbl in sorted(df_label['submodule'].unique()):
		inds = df_label['submodule']==lbl
		if np.sum(inds) > 10:
			smod_cnt += 1
			df_label.loc[inds, 'submodule'] = str(smod_cnt)
		else:
			df_label.loc[inds, 'submodule'] = ''

	return cg, df_label


if __name__ == '__main__':
	snrna_dir = '../snrna/'
	#snrna_file = os.path.join(snrna_dir, 'adata_U54_BA9_snrna_sennet.h5ad')
	# Update with cluster assignments from Jason
	snrna_file = os.path.join(snrna_dir, 'adata_U54_BA9_snrna_sennet_mapped_clusters.h5ad')

	# Revert to raw counts
	adata_sn = sc.read_h5ad(snrna_file)
	if hasattr(adata_sn, 'raw'):
		adata_sn = adata_sn.raw.to_adata()
	adata_sn.var.index = adata_sn.var['gene_name']

	adata_sn = preprocess_snrna(adata_sn)

	# Read module assignments and map to genes
	module_file = 'U54_BA9_splotch_lambdas_corr_hclust_labels_d1700.csv'
	df_modules = pd.read_csv(module_file, index_col=0)
	df_modules.dropna(axis=0, inplace=True)

	# Limit to set of genes in both tables (snRNA-seq, ST modules)
	genes_shared = np.intersect1d(adata_sn.var.index.values, df_modules['gene_name'].values)
	df_modules = df_modules.loc[df_modules['gene_name'].isin(genes_shared), :]
	adata_sn = adata_sn[:, genes_shared]

	mset = 'd1700'
	plot_dir = 'submodules'
	for m in sorted(df_modules[mset].unique()):
		if m==1:
			display_genes = ['PDGFRA', 'TIMP1', 'PECAM1', 'HEG1', 'MOG', 'MOBP', 
							 'PLP1', 'MBP', 'CNP', 'MAG', 'AIF1', 'TYROBP', 'CD68', 'C1QA']
		elif m==5:
			display_genes = ['CD34', 'DOCK9', 'FLT1', 'BCAM', 'ALDH1L1', 'SOX9', 'BMPR1B']
		elif m==13:
			display_genes = ['GLUL', 'TGFB2', 'LGALS3']
		else:
			display_genes = None
					   

		for obs_ctype in ['mapped_cluster']:
			if not os.path.exists(os.path.join(plot_dir, obs_ctype)):
				os.mkdir(os.path.join(plot_dir, obs_ctype))

			genes_m = df_modules['gene_name'][df_modules[mset]==m]

			if len(genes_m) > 10:
				g, df_label = plot_submodules(adata_sn, genes_m, obs_ctype=obs_ctype, display_genes=display_genes)
				
				g.ax_heatmap.set_title('Module %d (N=%d genes)' % (m, len(genes_m)), fontsize=22)
				plt.savefig(os.path.join(plot_dir, obs_ctype, 'd1700_m%d.png' % m), dpi=300, bbox_inches='tight')
				plt.close()

				df_label.to_csv(os.path.join(plot_dir, obs_ctype, 'd1700_m%d.csv' % m))
