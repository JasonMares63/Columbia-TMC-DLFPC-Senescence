import numpy as np
import pandas as pd
import scanpy as sc

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-l', '--lambda-file', type=str, required=True,
	help='Path to HDF5-formatted AnnData object with lambdas.')
parser.add_argument('-x', '--layer', type=str,
	help='Layer in AnnData object containing values on which to compute Pearson correlation.')
parser.add_argument('-d', '--dest-file', type=str, required=True,
	help='Path to save CSV-formatted DataFrame with gene-gene correlation matrix.')
args = parser.parse_args()

adata = sc.read_h5ad(args.lambda_file)

if args.layer is None:
        cap_val = np.percentile(adata.X, 99, axis=0)
        print('99th percentile:', cap_val, flush=True)
        lambdas_clipped = np.clip(adata.X, np.zeros_like(cap_val), cap_val)
        df_lambdas = pd.DataFrame(data=lambdas_clipped, columns=adata.var.index)
else:
        df_lambdas = pd.DataFrame(data=np.array(adata.layers[args.layer]), columns=adata.var.index)

corr_mat = df_lambdas.corr(method='pearson')
corr_mat.to_csv(args.dest_file)
