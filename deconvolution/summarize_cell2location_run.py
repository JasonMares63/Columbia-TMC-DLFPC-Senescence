import os
import glob
import pandas as pd
import scanpy as sc

deconv_results = []

for afile in glob.glob(os.path.join('.', '*', '*_cell2loc_deconv.h5ad')):
    adata = sc.read_h5ad(afile)
    
    '''
    # add 5% quantile, representing confident cell abundance, 'at least this amount is present'
    deconv_q5 = adata.obsm['q05_cell_abundance_w_sf']
    deconv_q5.columns = [x.split('abundance_w_sf_')[1] for x in deconv_q5.columns]
    deconv_results.append(deconv_q5)
    '''
    
    deconv = adata.obsm['means_cell_abundance_w_sf']
    deconv.columns = [x.split('abundance_w_sf_')[1] for x in deconv.columns]
    deconv_results.append(deconv)
    
deconv_results = pd.concat(deconv_results, axis=0)
deconv_results = deconv_results.fillna(0)
deconv_results.to_csv('mapped_clusters_comp_cell2location_markers_fullref.tsv', sep='\t')