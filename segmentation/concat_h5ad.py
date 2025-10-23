import pandas as pd
import numpy as np 
import scanpy as sc 
import glob
import os

h5ad_list = glob.glob("/sc/arion/projects/ji_lab/DATA/MERFISH/segmentation/NS-Atlas/baysor_0.7.1/*/*/qc/*.h5ad")

adata_list = []
for h5ad_file in h5ad_list:
    adata = sc.read_h5ad(h5ad_file)
    adata_list.append(adata)

import anndata as ad
adata_combined = ad.concat(adata_list, join='outer')
adata_combined.obs_names_make_unique()
adata_combined.write_h5ad("/sc/arion/projects/ji_lab/DATA/MERFISH/segmentation/NS-Atlas/baysor_0.7.1/ns-atlas.baysor_segmentation_filtered.anndata_object.h5ad")