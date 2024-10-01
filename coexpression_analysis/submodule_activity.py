import os
import glob
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import ttest_ind
from statsmodels.stats.multitest import fdrcorrection

def pval_to_stars(pval):
	if pval < 1e-4:
		return '****'
	elif pval < 1e-3:
		return '***'
	elif pval < 1e-2:
		return '**'
	elif pval < 0.05:
		return '*'
	else:
		return ''

def plot_module_heatmap(adata, module_genes, rows='aars', cols='Level 1', row_order=None, col_order=None,
						col_comparisons=None, group_by=None, vmin=None, vmax=None, ax=None):
	'''
	Parameters:
	----------
	adata: AnnData
		expression data to be visualized (should be scaled)
	module_genes: iterable of str
		names of genes in coexpression module of interest (must match naming in adata.var.index)
	rows: str
		annotations in adata.obs to separate by on rows of heatmap
	cols: str
		annotations in adata.obs to separate by on cols of heatmap
	row_order: iterable of str or None
		list of annotation categories, in order, to appear on heatmap rows (can be used to subset data)
	col_order: iterable of str or None
		list of annotation categories, in order, to appear on heatmpa cols (can be used to subset data)
	col_comparisons: iterable of tuple or None
		list of (col1, col2) comparisons across which to perform BH-adjusted (Welch's) t-test; denote on heatmap with *'s
	group_by: str
		annotations in adata.obs denoting independent groups of samples (e.g., Visium array names); defaults to all observations being independent
	vmin, vmax: float
		min/max values of color bar
	ax: Axes or None
		axes on which to plot heatmap, or None to instantiate new
		
	Returns:
	-------
	ax: Axes
	'''
	adata_sub = adata[:, module_genes]
	
	if row_order is None:
		row_order = adata.obs[rows].unique()
	if col_order is None:
		col_order = adata.obs[cols].unique()
		
	dat = dict([(c,[]) for c in col_order])
			
	# Calculate mean expression of module in each cell
	for c in col_order:
		adata_sub_c = adata_sub[adata_sub.obs[cols] == c]
		for r in row_order:
			adata_sub_r = adata_sub_c[adata_sub_c.obs[rows] == r]
			mean_expr = adata_sub_r.X.mean()
			dat[c].append(mean_expr)
	dat = pd.DataFrame(dat, index=row_order, dtype=np.float32)
	
	# Calculate significance of mean change between indicated column pairs (separately per row)
	sig = pd.DataFrame(data=np.ones_like(dat), index=dat.index, columns=dat.columns)
	for c1, c2 in col_comparisons:
		adata_sub_c1 = adata_sub[adata_sub.obs[cols] == c1]
		adata_sub_c2 = adata_sub[adata_sub.obs[cols] == c2]
		for r in row_order:
			adata_sub_r1 = adata_sub_c1[adata_sub_c1.obs[rows] == r]
			adata_sub_r2 = adata_sub_c2[adata_sub_c2.obs[rows] == r]
			
			means_1 = adata_sub_r1.X.mean(axis=1)
			means_2 = adata_sub_r2.X.mean(axis=1)
			# group averages (e.g., Visium arrays) treated as independent samples
			if group_by is not None:
				means_1 = [means_1[adata_sub_r1.obs[group_by]==g].mean() for g in adata_sub_r1.obs[group_by].unique()]
				means_2 = [means_2[adata_sub_r2.obs[group_by]==g].mean() for g in adata_sub_r2.obs[group_by].unique()]

			t, pval = ttest_ind(means_1, means_2, equal_var=False)
			sig.loc[r, c2] = pval
	
	# Perform BH FDR correction
	pvals_all = sig.values.flatten()
	pvals_all[pvals_all != 1.0] = fdrcorrection(pvals_all[pvals_all != 1.0], is_sorted=False)[1]
	p_adj = pd.DataFrame(pvals_all.reshape(sig.shape), index=sig.index, columns=sig.columns)
	p_stars = p_adj.applymap(pval_to_stars)
		
	if ax is None:
		fig, ax = plt.subplots(1)
	sns.heatmap(dat, ax=ax, cmap='bwr', vmin=vmin, vmax=vmax,
				annot=p_stars, fmt="",
				cbar_kws={'label':r'Mean scaled expression ($\overline{\lambda}$)'})
	
	# Add brackets over compared columns
	rendered = []
	for (c1, c2) in col_comparisons:
		if c1 in rendered or c2 in rendered:
			bh = -0.2
		else:
			bh = -0.1
		x1 = list(col_order).index(c1) + 0.5
		x2 = list(col_order).index(c2) + 0.5
		ymin, ymax = ax.get_ylim()
		ax.plot([x1, x1, x2, x2], [0, bh, bh, 0], c='k')
		ax.set_ylim(ymin, ymax+bh)

		rendered.append(c1)
		rendered.append(c2)

	xticklabels = ax.get_xticklabels()
	ax.set_xticklabels(xticklabels, fontsize=16)
	yticklabels = ax.get_yticklabels()
	ax.set_xticklabels(yticklabels, fontsize=16)
	
	return ax

if __name__ == '__main__':
	adata = sc.read_h5ad('../adata_U54_BA9_splotch_lambdas.h5ad')
	adata.obs['array'] = [x.split('/')[-2] for x in adata.obs.index]  # denote arrays in .obs

	# Use common gene names to match submodules
	gene_symbols = pd.read_csv('../gene_symbols.tsv', header=0, index_col=0, sep='\t')
	adata.var = adata.var.join(gene_symbols, how='left')
	adata = adata[:, adata.var.dropna().index].copy()  # drop genes without a common name
	adata.var.set_index('gene_name', inplace=True)

	# Remove effects of gene scale, so we can calculate mean submodule activities without certain genes dominating
	sc.pp.scale(adata)
	
	plot_dir = os.path.join('submodules', 'mapped_cluster_activities')
	if not os.path.exists(plot_dir):
		os.mkdir(plot_dir)

	for smfile in glob.glob(os.path.join('submodules', 'mapped_cluster', '*_m*.csv')):
		df = pd.read_csv(smfile, header=0, index_col=0)
		df = df.dropna().astype(int)

		m = int(smfile.split('_m')[-1].split('.csv')[0])

		for sm in df['submodule'].unique():
			sm_genes = df.index[df['submodule']==sm]

			ax = plot_module_heatmap(adata, sm_genes, rows='region', cols='Level 1',
				row_order=['Layer_1', 'Layer_2', 'Layer_3', 'Layer_4', 'Layer_5', 'Layer_6', 'White_matter'],
				col_order=['Young', 'Middle', 'Old'],
				col_comparisons=[('Young', 'Middle'), ('Young', 'Old')],
				group_by='array',
				vmin=-1, vmax=1
				)
			ax.set_title('Submodule %d.%d (n=%d genes)' % (m, sm, len(sm_genes)), fontsize=22)
			plt.savefig(os.path.join(plot_dir, 'submodule_%d-%d.png' % (m, sm)), dpi=300)
