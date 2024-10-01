import numpy as np
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path
from argparse import ArgumentParser
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster


parser = ArgumentParser()
parser.add_argument('-c', '--corr-mat', type=str, required=True,
	help='Path to CSV-formatted file containing gene-gene correlation matrix.')
parser.add_argument('-d', '--dest-file', type=str, required=True,
	help='Path to save linkage matrix output from hierarchical clustering.')
args = parser.parse_args()

df_corr = pd.read_csv(args.corr_mat, header=0, index_col=0)
# Pearson correlation can be NaN when one quantity has a constant value -- still means no correlation.
df_corr = df_corr.fillna(0)  


# At each iteration (starting from singleton clusters), merge the two clusters with
#   the smallest L1 distance between their centroids.
Z = linkage(df_corr.values, method='average', metric='cityblock')

np.save(Path(args.dest_file).with_suffix('.npy'), Z)

# Plot the full dendrogram so we can determine number of clusters
fig, ax = plt.subplots(1, figsize=(20,8))
dn1 = dendrogram(Z, 
	truncate_mode=None,
	leaf_rotation=90, leaf_font_size=6, ax=ax)
ax.set_xlabel('Gene index')
ax.set_ylabel('L1 distance')
plt.savefig(Path(args.dest_file).with_suffix('.jpg'), format='JPEG', dpi=300)


# Plot the correlation matrix with entries (genes) ordered according to hierarchical
#   cluster membership.
sns.clustermap(df_corr, method='average', metric='euclidean', 
		row_cluster=True, row_linkage=Z,
		col_cluster=True, col_linkage=Z,
		cmap='PRGn', 
		cbar_kws={'label':'Pearson correlation'},
		xticklabels=False, 
		yticklabels=False
	)
plt.savefig(Path(args.dest_file).with_suffix('.png'), format='PNG', dpi=300)