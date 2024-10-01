import os
import glob
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

from submodule_activity import plot_module_heatmap


adata = sc.read_h5ad('../adata_U54_BA9_splotch_lambdas.h5ad')
adata.obs['array'] = [x.split('/')[-2] for x in adata.obs.index]  # denote arrays in .obs

# Remove effects of gene scale, so we can calculate mean submodule activities without certain genes dominating
sc.pp.scale(adata)


# Read in external gene lists
df_external = pd.read_csv('SuppFig5_external_genelists.csv', header=0, skiprows=1)

output_dir = 'external_modules'
if not os.path.exists(output_dir):
	os.mkdir(output_dir)

for module in df_external.columns:
	genes_in = df_external[module].dropna()
	print(module, str(len(genes_in)), 'genes')
	genes_in = adata.var.index.intersection(genes_in).values
	print(len(genes_in), 'shared with Splotch genes')
	
	ax = plot_module_heatmap(adata, genes_in, rows='region', cols='Level 1',
		row_order=['Layer_1', 'Layer_2', 'Layer_3', 'Layer_4', 'Layer_5', 'Layer_6', 'White_matter'],
		col_order=['Young', 'Middle', 'Old'],
		col_comparisons=[('Young', 'Middle'), ('Young', 'Old')],
		group_by='array',
		vmin=-1, vmax=1
		)
	ax.set_title('%s (n=%d genes)' % (module, len(genes_in)))
	filename = module.replace(' ', '_')
	plt.savefig(os.path.join(output_dir, '%s.png' % filename), dpi=300)	
