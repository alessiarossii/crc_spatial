#!/usr/bin/env nextflow

process make_output {
    clusterOptions { "--job-name=nf-stAI-output-[${task.hash}]" }

    input:
        path(batches)

    output:
        path("spdata_imputed.h5ad")
    
    // Salva automaticamente i file nella cartella di output definita nei parametri NUOVO
    publishDir "${params.output_dir}/imputed", mode: 'copy'
    
    script:
	"""
#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import sys
import anndata as ad
import scanpy as sc
from scipy import sparse
import yaml

def parse_null(x):
    return None if x == 'null' else x

batch_key = parse_null('${params.batch_key}')
if batch_key is None:
    batch_key = 'single_batch_key'

batches = [ ${batches.collect{"'${it}'"}.join(', ')} ]

adata_list = []
pattern='spdata_imputed.h5ad'
for batch_dir in batches:
    for dirpath, dirnames, filenames in os.walk(batch_dir):
        if pattern in filenames:
            fpath = os.path.join(dirpath, pattern)
            print(f"Loading {fpath}")
            adata = sc.read_h5ad(fpath)
            adata_list.append(adata)

if not adata_list:
    raise ValueError(f"No files matching pattern '{pattern}' found in {batches}")
merged_adata = ad.concat(adata_list)
merged_adata.obs.index.name = None
merged_adata.obs_names_make_unique()
merged_adata.write_h5ad('spdata_imputed.h5ad', compression='gzip')

	"""

}